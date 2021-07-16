#include "common.h"
#include "spheres.h"

#include "study_case.h"
#include "case_tension_pointwise.h"
#include "case_tension.h"
#include "case_torsion.h"
#include "case_bending_displacement.h"
#include "case_bending_rotation.h"

#include "core/analysis_transient.h"
#include "core/solver_direct.h"
#include "core/solver_newton_raphson.h"
#include "core/timestep_selector_fixed_factor.h"
#include "core/action_matrix_solver.h"
#include "core/action_pipeline.h"
#include "core/triangulation_output.h"
#include "core/convergence_tester_residual.h"
#include "core/assembly_strategy_sequential.h"
#include "core/point_data_map.h"
#include "core/linear_soe_basic.h"
#include "core/analysis_one_step.h"

#include "core/scalar_quantity_time.h"
#include "core/scalar_quantity_const.h"
#include "core/table_output.h"
#include "core/table_entry_scalar.h"

#include "kernels/time_derivative.h"
#include "kernels/diffusion.h"
#include "kernels/functional_reaction.h"
#include "kernels/boundary_diffusion.h"

#include "pde/postprocessor.h"
#include "pde/quantity_variable_value.h"
#include "pde/quantity_scalar_material_object.h"
#include "pde/sensor_scalar.h"
#include "pde/field_interpolator.h"
#include "pde/builders.h"

#include "solid_mechanics/builders.h"
#include "solid_mechanics/material_plane_strain.h"
#include "solid_mechanics/material_tensor_isotropic.h"

#include <deal.II/base/function_lib.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/mapping_fe.h>

#include "linear_solver_sp.h"

// Main model itself
template <int dim>
Spheres<dim>::Spheres(dealii::ParameterHandler &param)
{
    // Read material properties
    param.enter_subsection("Physical constants");
    double E = param.get_double("E");
    double nu = param.get_double("nu");
    //bool isTensor = param.get_double("tensor");
    param.leave_subsection();

    // Read spheres centers
    param.enter_subsection("Geometry");
    // x1, x2
    double x1 = param.get_double("x1");
    double x2 = param.get_double("x2");

    dealii::Point<dim> O1;
    O1(0) = x1;
    
    dealii::Point<dim> O2;
    O2(0) = x2;

    TOOLS::cout::instance() << "x1 = " << x1 << std::endl;
    TOOLS::cout::instance() << "x2 = " << x2 << std::endl;

    const std::string mesh_file = param.get("mesh");
    const std::string stiffnessMode = param.get("stiffness");
    bool isRightActive = param.get_bool("right_active");

    TOOLS::cout::instance() << "Analysis case = " << stiffnessMode << std::endl;

    _triangulation = loadMesh(mesh_file);

    param.leave_subsection();

    // Linear solver - common for all problems
    unsigned int maxIter = 100000;
    double absTol = 1e-7;
    double relTol = 1e-7;
    auto linearSolver = std::make_shared<LinearSolverSp<VectorType, MatrixType>>(maxIter, absTol, relTol);

    // Plane elastic model
    _model = std::make_shared<PDE::ModelWithKernels<dim, VectorType>>();

    auto components = SOLID::createVectorVariables<dim, VectorType>(_model, "u");

    // Build kernels
    auto material = std::make_shared<SOLID::MaterialTensorIsotropic<dim>>(E, nu);
    SOLID::createTensorLinearKernels<dim, VectorType>(_model, material, components);

    // System assembly strategy - separate for each model
    using ModelType = std::remove_reference_t<decltype(*_model)>;
    auto assemblyStrategy
        = std::make_shared<CORE::AssemblyStrategySequential<dim, VectorType, MatrixType, ModelType>>(_model);

    // We create linear SOE explicitly
    auto linearSOE = std::make_shared<CORE::LinearSOEBasic<VectorType, MatrixType>>();

    // Nonlinear solver
    bool isNonlinear = false;
    std::shared_ptr<CORE::NonlinearSolver<dim, VectorType, MatrixType>> nonlinearSolver;
    if (isNonlinear) {
        // Newton solver
        double absTolNewton = 1e-8;
        unsigned int iterMax = 10;
        auto convergenceTester = std::make_shared<CORE::ConvergenceTesterResidual<VectorType, MatrixType>>(linearSOE, absTolNewton);
        nonlinearSolver = std::make_shared<CORE::SolverNewtonRaphson<dim, VectorType, MatrixType>>(iterMax, assemblyStrategy, linearSolver, convergenceTester);
    } else {
        // Linear one-step solver
        nonlinearSolver = std::make_shared<CORE::SolverDirect<dim, VectorType, MatrixType>>(assemblyStrategy, linearSolver);
    }
    /*
    std::shared_ptr<CORE::NonlinearSolver<dim, VectorType, MatrixType>> nonlinearSolver
        = std::make_shared<CORE::SolverDirect<dim, VectorType, MatrixType>>(assemblyStrategy, linearSolver);
    */
    
    // Create problem set
    _problemSet = std::make_shared<CORE::ProblemSet<dim, VectorType>>();
    _manager = std::make_shared<CORE::ActionsManager<dim, VectorType>>();

    // Approximation components for trias
    unsigned int degree = 1;

    dealii::FE_SimplexP<dim> feDummy(degree);
    dealii::FE_SimplexP<dim> fe(degree);
    auto quadCell = std::make_shared<dealii::QGaussSimplex<dim>>(degree + 1);
    auto quadFace = std::make_shared<dealii::QGaussSimplex<dim-1>>(degree + 1);

    dealii::FE_SimplexP<dim> fe_mapping(1);
    auto mapping = std::make_shared<dealii::MappingFE<dim>>(fe_mapping);

    // Create problems
    _problem = std::make_shared<CORE::Problem<dim, VectorType>>("Spheres", _model, quadCell, quadFace,
        _triangulation, _problemSet, std::move(fe), mapping);

    std::shared_ptr<StudyCase<dim, VectorType>> studyCase;
    if (stiffnessMode == "tension") {
        /*
        studyCase = std::make_shared<CaseTensionPointwise<dim, VectorType, MatrixType, ModelType>>(
            _problem, _model, linearSOE, O1, O2);
        */
        studyCase = std::make_shared<CaseTension<dim, VectorType, MatrixType, ModelType>>(
            _problem, _model, linearSOE, O1, O2, isRightActive);
    } else if (stiffnessMode == "torsion") {
        studyCase = std::make_shared<CaseTorsion<dim, VectorType, MatrixType, ModelType>>(
            _problem, _model, linearSOE, O1, O2, isRightActive);
    } else if (stiffnessMode == "bending_displacement") {
        studyCase = std::make_shared<CaseBendingDisplacement<dim, VectorType, MatrixType, ModelType>>(
            _problem, _model, linearSOE, O1, O2, isRightActive);
    } else if (stiffnessMode == "bending_rotation") {
        studyCase = std::make_shared<CaseBendingRotation<dim, VectorType, MatrixType, ModelType>>(
            _problem, _model, linearSOE, O1, O2, isRightActive);
    }
    studyCase->imposeBoundaryConditions(_model);

    // Apply default kinematic EBCs
    //setBoundaryCondtions(_model, *_triangulation, O1, O2, stiffnessMode);

    // 1. Solve
    _manager->appendAction( std::make_shared<CORE::ActionMatrixSolver<dim, VectorType, MatrixType>>(_problem, nonlinearSolver, linearSOE) );

    // Output manager
    auto outputManager = std::make_shared<CORE::OutputManager>();
    
    param.enter_subsection("Output");
    const std::string savePathVtk = param.get("path_vtk");
    const std::string savePathLog = param.get("path_log");
    bool doSaveVtk = param.get_bool("save_vtk");
    bool doSaveLog = param.get_bool("save_log");
    param.leave_subsection();

    // Output objects
    if (doSaveVtk) {
        auto output = std::make_shared<CORE::TriangulationOutput<dim, VectorType>>(_triangulation, savePathVtk);
        auto postprocessor = std::make_shared<PDE::Postprocessor<dim, VectorType>>(
            SOLID::buildMechanicsPostprocessor<dim, VectorType>(_model, _problem, components));
        output->addPostprocessor(postprocessor);

        // Add output
        outputManager->addOutput(output);
    }

    if (doSaveLog) {
        // Table output for the primal model
        auto table = std::make_shared<CORE::TableOutput>(savePathLog);

        // Add table output
        outputManager->addOutput(table);

        auto scalarTime = std::make_shared<CORE::ScalarQuantityTime>(_model);
        table->addEntry(std::make_shared<CORE::TableEntryScalar<double>>(scalarTime));

        auto scalarLoadMagnitude = std::make_shared<CORE::ScalarQuantityConst>("load_magnitude", studyCase->loadMagnitude());
        table->addEntry(std::make_shared<CORE::TableEntryScalar<double>>(scalarLoadMagnitude));

        studyCase->addReactions(table);
    }

    // Timestep selector settings
    double timestepGrowthFactor = 1.0;
    auto timestepSelector = std::make_shared<CORE::TimestepSelectorFixedFactor>(timestepGrowthFactor);

    param.enter_subsection("Time integration");
    
    param.leave_subsection();

    // Analysis settings
    _analysis = std::make_shared<CORE::AnalysisOneStep<dim, VectorType>>(outputManager);
}

template <int dim>
void Spheres<dim>::run()
{
    _analysis->analyze(_manager);
}
/*
template <int dim>
std::shared_ptr<Triangulation<dim>> Spheres<dim>::loadMesh(const std::string& filename)
{
    auto triangulation = createTriaPtr<dim>();

    dealii::GridIn<dim> gridin;
    gridin.attach_triangulation(*triangulation);
    std::ifstream f(filename);
    gridin.read_msh(f);

    // Second mesh
    std::string filename2 = "/home/kapusta/work/ai-2particles/spheres2-v1-part2.msh";
    auto triangulation2 = createTriaPtr<dim>();
    dealii::GridIn<dim> gridin2;
    gridin2.attach_triangulation(*triangulation2);
    std::ifstream f2(filename2);
    gridin2.read_msh(f2);

    //dealii::GridIn<dim>(*triangulation).read(filename);
    const auto refere_cell_types = triangulation->get_reference_cells();

    AssertDimension(refere_cell_types.size(), 1);


    auto triangulationRes = createTriaPtr<dim>();

    dealii::GridGenerator::merge_triangulations(*triangulation, *triangulation2, *triangulationRes);

    // Output grid for debug

    std::ofstream out1("grid1.svg");
    dealii::GridOut grid_out1;
    grid_out1.write_svg(*triangulation, out1);

    std::ofstream out2("grid2.svg");
    dealii::GridOut grid_out2;
    grid_out2.write_svg(*triangulation2, out2);

    std::ofstream outRes("gridRes.svg");
    dealii::GridOut grid_outRes;
    grid_outRes.write_svg(*triangulationRes, outRes);

    throw std::runtime_error("STOP after Grid");

    return triangulationRes;
}
*/

template <int dim>
std::shared_ptr<Triangulation<dim>> Spheres<dim>::loadMesh(const std::string& filename)
{
    auto triangulation = createTriaPtr<dim>();

    dealii::GridIn<dim> gridin;
    gridin.attach_triangulation(*triangulation);
    std::ifstream f(filename);
    gridin.read_msh(f);
/*
    std::ofstream out("grid.svg");
    dealii::GridOut grid_out;
    grid_out.write_svg(*triangulation, out);

    throw std::runtime_error("STOP after Grid");
*/
    return triangulation;
}

/*
template <int dim>
void Spheres<dim>::setBoundaryCondtions(std::shared_ptr<CORE::Model<dim, VectorType>> model, Triangulation<dim>& triangulation, 
    const dealii::Point<dim>& O1, const dealii::Point<dim>& O2, const std::string& mode)
{
    (void)triangulation;

    if (mode == "tension") {
        model->addVertexConstraint(0, O1);
        model->addVertexConstraint(1, O1);
        if (dim == 3) {
            model->addVertexConstraint(2, O1);
        }

        model->addVertexConstraint(0, O2, 1);
        model->addVertexConstraint(1, O2);
        if (dim == 3) {
            model->addVertexConstraint(2, O2);
        }
    }

    if (mode == "torsion") {

        double posTol = 1e-6;

        std::map<dealii::Point<dim>, std::vector<unsigned int>> mapLeft;
        std::map<dealii::Point<dim>, std::vector<unsigned int>> mapRight;

        for (const auto &tria_cell : triangulation.active_cell_iterators())
        {
            if (tria_cell->is_locally_owned())
            {
                auto cellSrc = typename dealii::DoFHandler<dim>::cell_iterator(*tria_cell, 
                    &_problem->getApproximation().getDoFHandler());

                for (unsigned int v = 0; v < cellSrc->n_vertices(); ++v)
                {
                    auto& point = tria_cell->vertex(v);

                    if (std::abs(point(0) - O1(0)) < posTol) {
                        auto idY = cellSrc->vertex_dof_index(v, 1);
                        auto idZ = cellSrc->vertex_dof_index(v, 2);

                        mapLeft[point] = {idY, idZ}
                    } else if (std::abs(point(0) - O2(0)) < posTol) {
                        auto idY = cellSrc->vertex_dof_index(v, 1);
                        auto idZ = cellSrc->vertex_dof_index(v, 2);

                        mapRight[point] = {idY, idZ}
                    }
                }
            }
        }
    }
}
*/

template class Spheres<2>;
template class Spheres<3>;