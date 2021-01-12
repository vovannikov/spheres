#include "bone.h"
#include "components.h"
#include "common.h"

#include "core/analysis_transient.h"
#include "core/solver_direct.h"
#include "core/solver_newton_raphson.h"
#include "core/timestep_selector_fixed_factor.h"
#include "core/action_matrix_solver.h"
#include "core/action_pipeline.h"
#include "core/triangulation_output.h"
#include "core/convergence_tester_residual.h"
#include "core/linear_soe.h"
#include "core/assembly_strategy_sequential.h"
#include "core/point_data_map.h"

#include "kernels/time_derivative.h"
#include "kernels/diffusion.h"
#include "kernels/functional_reaction.h"
#include "kernels/boundary_diffusion.h"

#include "pde/material_functor_wrapper.h"
#include "pde/postprocessor.h"
#include "pde/quantity_variable_value.h"
#include "pde/quantity_scalar_material_object.h"
#include "pde/sensor_scalar.h"
#include "pde/field_interpolator.h"
#include "pde/builders.h"

#include <deal.II/base/function_lib.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/simplex/fe_lib.h>
#include <deal.II/simplex/grid_generator.h>
#include <deal.II/simplex/quadrature_lib.h>

#include <deal.II/fe/mapping_fe.h>

// Main model itself
template <int dim>
Bone<dim>::Bone(dealii::ParameterHandler &param)
{
    param.enter_subsection("Physical constants");
    // Parameters
    // DM1
    double D_mM1 = param.get_double("D_mM1");
    double K_mM1_D = param.get_double("K_mM1_D");

    // DM2
    double D_mM2 = param.get_double("D_mM2");
    double K_mM2_D = param.get_double("K_mM2_D");

    // CM1
    double C_TNFM1 = param.get_double("C_TNFM1");
    double K_TNFM1_C = param.get_double("K_TNFM1_C");
    double C_mM1 = param.get_double("C_mM1");
    double K_mM1_C = param.get_double("K_mM1_C");

    // AM1
    double A_TNFM1 = param.get_double("A_TNFM1");
    double A_sM1 = param.get_double("A_sM1");
    double K_TNFM1_A = param.get_double("K_TNFM1_A");
    double K_sM1_A = param.get_double("K_sM1_A");

    // AM2
    double A_IL4M2 = param.get_double("A_IL4M2");
    double A_sM2 = param.get_double("A_sM2");
    double K_IL4M2_A = param.get_double("K_IL4M2_A");
    double K_sM2_A = param.get_double("K_sM2_A");

    // F5
    double I_IL4M1 = param.get_double("I_IL4M1");
    double I_TGFM1 = param.get_double("I_TGFM1");
    double K_IL4M1_I = param.get_double("K_IL4M1_I");
    double K_TGFM1_I = param.get_double("K_TGFM1_I");

    // YM1
    double Y_M1_h = param.get_double("Y_M1_h");
    double Y_M1_L = param.get_double("Y_M1_L");
    double K_TNFM1_Y = param.get_double("K_TNFM1_Y");

    // YM2
    double Y_M2_h = param.get_double("Y_M2_h");

    // ETNF
    double E_TNF_h = param.get_double("E_TNF_h");
    double K_TNF = param.get_double("K_TNF");

    // EIL4
    double E_IL4_h = param.get_double("E_IL4_h");
    double K_IL4 = param.get_double("K_IL4");

    // DTNF
    double D_TNF_h = param.get_double("D_TNF_h");

    // DIL4
    double D_IL4_h = param.get_double("D_IL4_h");

    // HTNF
    double d_TNF = param.get_double("d_TNF");

    // HIL4
    double d_IL4 = param.get_double("d_IL4");

    // Dd
    double D_d = param.get_double("D_d");

    // dd
    double d_d = param.get_double("d_d");

    param.leave_subsection();


    param.enter_subsection("Geometry");

    const std::string mesh_file = param.get("mesh");

    //_triangulation = createMesh();
    _triangulation = loadMesh(mesh_file);

    param.leave_subsection();

    // Linear solver - common for all problems
    unsigned int maxIter = 1000;
    double absTol = 1e-8;
    double relTol = 1e-8;
    auto linearSolver = std::make_shared<LinearSolverType<VectorType, MatrixType>>(maxIter, absTol, relTol);

    auto modelDummy = std::make_shared<PDE::ModelWithVariables<dim, VectorType>>();
    _model = std::make_shared<PDE::ModelWithKernels<dim, VectorType>>();

    // Dummy variable for mf
    auto mVar = modelDummy->getSystem().createVariable("mf", false);

    // Variables - 5
    auto cm1Var = _model->getSystem().createVariable("cm1", true);
    auto cm2Var = _model->getSystem().createVariable("cm2", true);
    auto gTNFVar = _model->getSystem().createVariable("gTNF", true);
    auto gIL4Var = _model->getSystem().createVariable("gIL4", true);
    auto mdVar = _model->getSystem().createVariable("md", true);

    _model->getSystem().appendExternalVariable(mVar);

    // Components
    auto cmpDM1 = std::make_shared<Di<dim>>(mVar, D_mM1, K_mM1_D);
    auto cmpDM2 = std::make_shared<Di<dim>>(mVar, D_mM2, K_mM2_D);
    auto cmpCM1 = std::make_shared<CM1<dim>>(gTNFVar, mVar, cm1Var, C_TNFM1, K_TNFM1_C, C_mM1, K_mM1_C);
    auto cmpAM1 = std::make_shared<Ai<dim>>(mVar, cm1Var, cm2Var, gTNFVar, A_TNFM1, A_sM1, K_TNFM1_A, K_sM1_A);
    auto cmpAM2 = std::make_shared<Ai<dim>>(mVar, cm1Var, cm2Var, gIL4Var, A_IL4M2, A_sM2, K_IL4M2_A, K_sM2_A);
    auto cmpF5 = std::make_shared<F5<dim>>(gIL4Var, gTNFVar, I_IL4M1, I_TGFM1, K_IL4M1_I, K_TGFM1_I); //gTGFVar??
    auto cmpYM1 = std::make_shared<YM1<dim>>(gTNFVar, Y_M1_h, Y_M1_L, K_TNFM1_Y);
    auto cmpYM2 = std::make_shared<Dconst<dim>>(Y_M2_h);
    auto cmpETNF = std::make_shared<ETNF<dim>>(mdVar, gTNFVar, E_TNF_h, K_TNF);
    auto cmpEIL4 = std::make_shared<EIL4<dim>>(gIL4Var, E_IL4_h, K_IL4);
    auto cmpDTNF = std::make_shared<Dconst<dim>>(D_TNF_h);
    auto cmpDIL4 = std::make_shared<Dconst<dim>>(D_IL4_h);
    auto cmpHTNF = std::make_shared<H<dim>>(cm1Var, cm2Var, d_TNF);
    auto cmpHIL4 = std::make_shared<H<dim>>(cm1Var, cm2Var, d_IL4);
    auto cmpDd = std::make_shared<Dconst<dim>>(D_d);
    auto cmpdd = std::make_shared<Dconst<dim>>(d_d);

    // 1.1 Time derivative - cm1
    auto cm1KrnlTimeDer = std::make_shared<PDE::KernelTimeDerivative<dim>>(cm1Var, _model);
    _model->getSystem().addKernel(cm1KrnlTimeDer);

    // 1.2 Diffusion for cm1 - cm1
    auto cm1KrnlDiffCm1 = std::make_shared<PDE::KernelDiffusion<double, dim>>(cm1Var, cmpDM1);
    _model->getSystem().addKernel(cm1KrnlDiffCm1);

    // 1.3 Diffusion for gTNF - cm1
    auto cm1KrnlDiffgTNF = std::make_shared<PDE::KernelDiffusion<double, dim>>(cm1Var, gTNFVar, cmpCM1);
    _model->getSystem().addKernel(cm1KrnlDiffgTNF);

    // 1.4 Reaction AM1 - cm1
    auto cm1KrnlReacAM1 = std::make_shared<PDE::KernelFunctionalReaction<double, dim>>(cm1Var, cmpAM1);
    _model->getSystem().addKernel(cm1KrnlReacAM1);
    cm1KrnlReacAM1->coupleWith(cm2Var);
    cm1KrnlReacAM1->coupleWith(gTNFVar);

    // 1.5 Reaction F5 - cm1
    auto cm1KrnlReacF5 = std::make_shared<PDE::KernelFunctionalReaction<double, dim>>(cm1Var, cmpF5);
    _model->getSystem().addKernel(cm1KrnlReacF5);
    cm1KrnlReacF5->coupleWith(gIL4Var);
    cm1KrnlReacF5->coupleWith(gTNFVar);

    // 1.6 Reaction YM1 - cm1
    auto cm1KrnlReacYM1 = std::make_shared<PDE::KernelFunctionalReaction<double, dim>>(cm1Var, cmpYM1);
    _model->getSystem().addKernel(cm1KrnlReacYM1);
    cm1KrnlReacYM1->coupleWith(gTNFVar);

    //---------------------------------------------------------

    // 2.1 Time derivative - cm2
    auto cm2KrnlTimeDer = std::make_shared<PDE::KernelTimeDerivative<dim>>(cm2Var, _model);
    _model->getSystem().addKernel(cm2KrnlTimeDer);

    // 2.2 Diffusion for cm2 - cm2
    auto cm2KrnlDiffCm2 = std::make_shared<PDE::KernelDiffusion<double, dim>>(cm2Var, cmpDM2);
    _model->getSystem().addKernel(cm2KrnlDiffCm2);

    // 2.3 Reaction AM2 - cm2
    auto cm2KrnlReacAM2 = std::make_shared<PDE::KernelFunctionalReaction<double, dim>>(cm2Var, cmpAM2);
    _model->getSystem().addKernel(cm2KrnlReacAM2);
    cm2KrnlReacAM2->coupleWith(cm1Var);
    cm2KrnlReacAM2->coupleWith(gIL4Var);

    // 2.4 Reaction F5 from cm1 - cm2
    auto cmpF5minus = std::make_shared<PDE::MaterialFunctorWrapper<double, dim>>(cmpF5, [](double val) { return -val; });
    auto cm2KrnlReacF5 = std::make_shared<PDE::KernelFunctionalReaction<double, dim>>(cm2Var, cm1Var, cmpF5minus);
    _model->getSystem().addKernel(cm2KrnlReacF5);
    cm2KrnlReacF5->coupleWith(gIL4Var);
    cm2KrnlReacF5->coupleWith(gTNFVar);

    // 2.5 Reaction YM2 - cm2
    auto cm2KrnlReacYM2 = std::make_shared<PDE::KernelFunctionalReaction<double, dim>>(cm2Var, cmpYM2);
    _model->getSystem().addKernel(cm2KrnlReacYM2);

    //---------------------------------------------------------

    // 3.1 Time derivative - gTNF
    auto gTNFKrnlTimeDer = std::make_shared<PDE::KernelTimeDerivative<dim>>(gTNFVar, _model);
    _model->getSystem().addKernel(gTNFKrnlTimeDer);

    // 3.2 Diffusion DTNF - gTNF
    auto gTNFKrnlDiffgTNF = std::make_shared<PDE::KernelDiffusion<double, dim>>(gTNFVar, cmpDTNF);
    _model->getSystem().addKernel(gTNFKrnlDiffgTNF);

    // 3.3 Reaction ETNF - gTNF
    auto gTNFKrnlReacETNF = std::make_shared<PDE::KernelFunctionalReaction<double, dim>>(gTNFVar, cm1Var, cmpETNF);
    _model->getSystem().addKernel(gTNFKrnlReacETNF);
    gTNFKrnlReacETNF->coupleWith(mdVar);

    // 3.4 Reaction H from cm1 + cm2 - gTNF
    auto gTNFKrnlReacHTNF = std::make_shared<PDE::KernelFunctionalReaction<double, dim>>(gTNFVar, cmpHTNF);
    _model->getSystem().addKernel(gTNFKrnlReacHTNF);
    gTNFKrnlReacHTNF->coupleWith(cm1Var);
    gTNFKrnlReacHTNF->coupleWith(cm2Var);

    //---------------------------------------------------------

    // 4.1 Time derivative - gIL4
    auto gIL4KrnlTimeDer = std::make_shared<PDE::KernelTimeDerivative<dim>>(gIL4Var, _model);
    _model->getSystem().addKernel(gIL4KrnlTimeDer);

    // 4.2 Diffusion DIL4 - gIL4
    auto gIL4KrnlDiffgIL4 = std::make_shared<PDE::KernelDiffusion<double, dim>>(gIL4Var, cmpDIL4);
    _model->getSystem().addKernel(gIL4KrnlDiffgIL4);

    // 4.3 Reaction EIL4 - gIL4
    auto gIL4KrnlReacEIL4 = std::make_shared<PDE::KernelFunctionalReaction<double, dim>>(gIL4Var, cm2Var, cmpEIL4);
    _model->getSystem().addKernel(gIL4KrnlReacEIL4);

    // 4.4 Reaction H from cm1 + cm2 - gIL4
    auto gIL4KrnlReacHIL4 = std::make_shared<PDE::KernelFunctionalReaction<double, dim>>(gIL4Var, cmpHIL4);
    _model->getSystem().addKernel(gIL4KrnlReacHIL4);
    gIL4KrnlReacHIL4->coupleWith(cm1Var);
    gIL4KrnlReacHIL4->coupleWith(cm2Var);

    //---------------------------------------------------------

    // 5.1 Time derivative - md
    auto mdKrnlTimeDer = std::make_shared<PDE::KernelTimeDerivative<dim>>(mdVar, _model);
    _model->getSystem().addKernel(mdKrnlTimeDer);

    // 5.2 Diffusion Dd - md
    auto mdKrnlDiffMd = std::make_shared<PDE::KernelDiffusion<double, dim>>(mdVar, cmpDd);
    _model->getSystem().addKernel(mdKrnlDiffMd);

    // 5.3 Reaction dd for cm1 - md
    auto mdKrnlReacCm1 = std::make_shared<PDE::KernelFunctionalReaction<double, dim>>(mdVar, cm1Var, cmpdd);
    _model->getSystem().addKernel(mdKrnlReacCm1);

    // 5.4 Reaction dd for cm2 - md
    auto mdKrnlReacCm2 = std::make_shared<PDE::KernelFunctionalReaction<double, dim>>(mdVar, cm2Var, cmpdd);
    _model->getSystem().addKernel(mdKrnlReacCm2);



    // Apply periodic boundary conditions
    setBoundaryCondtions(_model, *_triangulation);

    // System assembly strategy - separate for each model
    auto assemblyStrategy
        = std::make_shared<CORE::AssemblyStrategySequential<dim, VectorType, MatrixType, std::remove_reference_t<decltype(*_model)>>>(_model);

    // We create linear SOE explicitly
    auto linearSOE = std::make_shared<CORE::LinearSOE<VectorType, MatrixType>>();

    // Nonlinear solver
    std::shared_ptr<CORE::NonlinearSolver<dim, VectorType, MatrixType>> nonlinearSolver;
    if (_isNonlinear) {
        // Newton solver
        double absTolNewton = 1e-8;
        unsigned int iterMax = 10;
        auto convergenceTester = std::make_shared<CORE::ConvergenceTesterResidual<VectorType, MatrixType>>(linearSOE, absTolNewton);
        nonlinearSolver = std::make_shared<CORE::SolverNewtonRaphson<dim, VectorType, MatrixType>>(iterMax, assemblyStrategy, linearSolver, convergenceTester);
    } else {
        // Linear one-step solver
        nonlinearSolver = std::make_shared<CORE::SolverDirect<dim, VectorType, MatrixType>>(assemblyStrategy, linearSolver);
    }

    // Create problem set
    _problemSet = std::make_shared<CORE::ProblemSet<dim, VectorType>>();
    _manager = std::make_shared<CORE::ActionsManager<dim, VectorType>>();

    // Approximation components for trias
    unsigned int degree = 1;

    dealii::Simplex::FE_P<dim> feDummy(degree);
    dealii::Simplex::FE_P<dim> fe(degree);
    auto quadCell = std::make_shared<dealii::Simplex::QGauss<dim>>(degree + 1);
    auto quadFace = std::make_shared<dealii::Simplex::QGauss<dim-1>>(degree + 1);

    dealii::Simplex::FE_P<dim> fe_mapping(1);
    auto mapping = std::make_shared<dealii::MappingFE<dim>>(fe_mapping);

    // Create problems
    auto problemDummy = std::make_shared<CORE::Problem<dim, VectorType>>("Dummy", modelDummy, quadCell, quadFace,
        _triangulation, _problemSet, std::move(feDummy), mapping);
    _problem = std::make_shared<CORE::Problem<dim, VectorType>>("Bone", _model, quadCell, quadFace,
        _triangulation, _problemSet, std::move(fe), mapping);

    // Initial conditions for the system
    param.enter_subsection("Initial conditions");

    double initMf = param.get_double("mf");

    // 1. Apply initial conditions
    auto initialMf = std::make_shared<dealii::Functions::ConstantFunction<dim>>(initMf, 1);
    _fieldInterpolator = std::make_shared<PDE::FieldInterpolator<dim, VectorType>>();
    _fieldInterpolator->applyInterpolation(initialMf);
    _manager->appendOpenTimeAction( std::make_shared<CORE::ActionPipeline<dim, VectorType>>(problemDummy, _fieldInterpolator) );

    double initCm1 = param.get_double("cm1");
    double initCm2 = param.get_double("cm2");
    double initGTNF = param.get_double("gTNF");
    double initGIL4 = param.get_double("gIL4");
    double initMd = param.get_double("md");
    
    param.leave_subsection();

    auto initialValues = std::make_shared<dealii::Functions::ConstantFunction<dim>>(
        std::vector<double>{initCm1, initCm2, initGTNF, initGIL4, initMd});
    _problem->setInitialConditions(initialValues);

    // 2. Solve
    _manager->appendAction( std::make_shared<CORE::ActionMatrixSolver<dim, VectorType, MatrixType>>(_problem, nonlinearSolver, linearSOE) );

    // Output manager
    auto outputManager = std::make_shared<CORE::OutputManager>();
    
    param.enter_subsection("Output");
    const std::string savePath = param.get("save_path");
    param.leave_subsection();

    if (!savePath.empty()) {
        // Output objects
        auto output = std::make_shared<CORE::TriangulationOutput<dim, VectorType>>(_triangulation, savePath, nullptr, mapping);

        // Output mf
        // Build dummy postprocessor
        auto postprocessorDummy = std::make_shared<PDE::Postprocessor<dim, VectorType>>(problemDummy, modelDummy);
        postprocessorDummy->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityVariableValue<dim>>(mVar)));
        output->addPostprocessor(postprocessorDummy);


        // Build main postprocessor
        auto postprocessor = std::make_shared<PDE::Postprocessor<dim, VectorType>>(_problem, _model);

        // Add variables to output
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityVariableValue<dim>>(cm1Var)));
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityVariableValue<dim>>(cm2Var)));
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityVariableValue<dim>>(gTNFVar)));
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityVariableValue<dim>>(gIL4Var)));
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityVariableValue<dim>>(mdVar)));

        // Add components to output
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityScalarMaterialObject<dim>>("TERM_DM1", cmpDM1)));
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityScalarMaterialObject<dim>>("TERM_DM2", cmpDM2)));
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityScalarMaterialObject<dim>>("TERM_CM1", cmpCM1)));
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityScalarMaterialObject<dim>>("TERM_AM1", cmpAM1)));
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityScalarMaterialObject<dim>>("TERM_AM2", cmpAM2)));
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityScalarMaterialObject<dim>>("TERM_F5", cmpF5)));
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityScalarMaterialObject<dim>>("TERM_YM1", cmpYM1)));
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityScalarMaterialObject<dim>>("TERM_YM2", cmpYM2)));
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityScalarMaterialObject<dim>>("TERM_ETNF", cmpETNF)));
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityScalarMaterialObject<dim>>("TERM_EIL4", cmpEIL4)));
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityScalarMaterialObject<dim>>("TERM_DTNF", cmpDTNF)));
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityScalarMaterialObject<dim>>("TERM_DIL4", cmpDIL4)));
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityScalarMaterialObject<dim>>("TERM_HTNF", cmpHTNF)));
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityScalarMaterialObject<dim>>("TERM_HIL4", cmpHIL4)));
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityScalarMaterialObject<dim>>("TERM_Dd", cmpDd)));
        postprocessor->addSensor(std::make_shared<PDE::SensorScalar<dim>>(std::make_shared<PDE::QuantityScalarMaterialObject<dim>>("TERM_dd", cmpdd)));

        // Add postprocessor
        output->addPostprocessor(postprocessor);

        // Add output
        outputManager->addOutput(output);
    }

    // Timestep selector settings
    double timestepGrowthFactor = 1.0;
    auto timestepSelector = std::make_shared<CORE::TimestepSelectorFixedFactor>(timestepGrowthFactor);

    param.enter_subsection("Time integration");

    _timeStart = param.get_double("time_start");
    _timeEnd = param.get_double("time_end");
    _timeStepInit = param.get_double("time_step_init");
    _timeStepMin = param.get_double("time_step_min");
    _timeStepMax = param.get_double("time_step_max");
    
    param.leave_subsection();

    // Analysis settings
    _analysis = std::make_shared<CORE::AnalysisTransient<dim, VectorType>>(outputManager, timestepSelector);
}

template <int dim>
void Bone<dim>::run()
{
    _analysis->analyze(_manager, _timeStart, _timeEnd, _timeStepInit, _timeStepMin, _timeStepMax);
}

template <int dim>
std::shared_ptr<Triangulation<dim>> Bone<dim>::loadMesh(const std::string& filename)
{
    auto triangulation = createTriaPtr<dim>();

    dealii::GridIn<2> gridin;
    gridin.attach_triangulation(*triangulation);
    std::ifstream f(filename);
    gridin.read_msh(f);

    //dealii::GridIn<2>(*triangulation).read(filename);
    const auto refere_cell_types = triangulation->get_reference_cell_types();

    AssertDimension(refere_cell_types.size(), 1);

    // Output grid for debug
    /*
    std::ofstream out("grid.svg");
    dealii::GridOut grid_out;
    grid_out.write_svg(*triangulation, out);

    throw std::runtime_error("STOP after Grid");
    */

    // Source of cm1
    for (const auto &cell : triangulation->active_cell_iterators()) {
        
        for (const auto &face : cell->face_iterators()) {

            const dealii::Point<dim> face_center = face->center();

            if (face->at_boundary() && std::fabs(face_center(0)) > 1e-8 && std::fabs(face_center(1)) > 1e-8) {
                face->set_boundary_id(1);
            }
        }
    }

    return triangulation;
}

template <int dim>
std::shared_ptr<Triangulation<dim>> Bone<dim>::createMesh()
{
    // Part 1
    std::vector<unsigned int> subdivisions1{23, 27};//{11, 13}
    const dealii::Point< dim > bottom_left1(0, 0);
    const dealii::Point< dim > top_right1(1.14, 1.35);

    auto triangulation1 = createTriaPtr<dim>();

    dealii::GridGenerator::subdivided_hyper_rectangle(*triangulation1,
        subdivisions1, bottom_left1, top_right1);

    // Part 2
    std::vector<unsigned int> subdivisions2{17, 10};//{8, 5}
    const dealii::Point< dim > bottom_left2(1.14, 0);
    const dealii::Point< dim > top_right2(1.14 + 0.85, 0.5);

    auto triangulation2 = createTriaPtr<dim>();

    dealii::GridGenerator::subdivided_hyper_rectangle(*triangulation2,
        subdivisions2, bottom_left2, top_right2);


    // Part 3
    std::vector<dealii::Point<dim>> vertices;
    vertices.emplace_back(1.14 + 0.85, 0);
    vertices.emplace_back(1.14 + 0.85 + 1.43, 0);
    vertices.emplace_back(1.14 + 0.85, 4.35);
    
    auto triangulation3 = createTriaPtr<dim>();

    dealii::GridGenerator::simplex(*triangulation3, vertices);

    // Merged tria
    auto min_line_length = [](const Triangulation<dim> &tria) -> double
    {
        double length = std::numeric_limits<double>::max();
        for (const auto &cell : tria.active_cell_iterators())
            for (unsigned int n = 0; n < dealii::GeometryInfo<dim>::lines_per_cell; ++n)
                length = std::min(length, (cell->line(n)->vertex(0) - cell->line(n)->vertex(1)).norm());
        return length;
    };
    const double tolerance = std::min(min_line_length(*triangulation1),
                                    min_line_length(*triangulation2)) / 2.0;

    auto triangulation = createTriaPtr<dim>();

    dealii::GridGenerator::merge_triangulations(*triangulation1, *triangulation2,
        *triangulation, tolerance);

    // Output grid for debug
    std::ofstream out("grid.svg");
    dealii::GridOut grid_out;
    grid_out.write_svg(*triangulation, out);

    throw std::runtime_error("STOP after Grid");

    /* Mark boundaries
                  3
            ------------
            |          |
          0 |          | 1
            |          |
            ------------ x
                  2
    */
    /*
    for (const auto &cell : triangulation->active_cell_iterators()) {
        for (unsigned int face = 0; face < dealii::GeometryInfo<dim>::faces_per_cell; ++face) {

            if (std::fabs(cell->face(face)->center()(0) - 0) < 1e-8) {
                // left face
                cell->face(face)->set_boundary_id(0);
            } else if (std::fabs(cell->face(face)->center()(0) - length) < 1e-8) {
                // right face
                cell->face(face)->set_boundary_id(1);
            } else if (std::fabs(cell->face(face)->center()(1) - 0) < 1e-8) {
                // bottom face
                cell->face(face)->set_boundary_id(2);
            } else if (std::fabs(cell->face(face)->center()(1) - length) < 1e-8) {
                // top face
                cell->face(face)->set_boundary_id(3);
            }
        }
    }
    */

    return triangulation;
}

template <int dim>
void Bone<dim>::setBoundaryCondtions(std::shared_ptr<CORE::Model<dim, VectorType>> model, Triangulation<dim>& triangulation)
{
    (void)triangulation;
    auto bcVals = std::make_shared<dealii::Functions::ConstantFunction<dim>>(std::vector<double>{1, 0, 0, 0, 0});
    model->addDirichletBoundary(1, dealii::ComponentMask(std::vector<bool>{true, false, false, false, false}), bcVals);
}

template class Bone<2>;
//template class Bone<3>;