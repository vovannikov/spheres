#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>

#include "pde/functional_object.h"
#include "pde/material_object.h"

#include "common.h"

// Forward declarations
namespace PDE {
    template <int dim, typename VectorType>
    class ModelWithKernels;

    template <int dim, typename VectorType>
    class FieldInterpolator;

    template <int dim>
    class VariableValues;

    template <int dim>
    class Variable;
}

namespace CORE {
    template <int dim, typename VectorType>
    class AnalysisOneStep;

    template <int dim, typename VectorType>
    class ProblemSet;

    template <int dim, typename VectorType>
    class Problem;

    template <int dim, typename VectorType>
    class ActionsManager;

    template <int dim, typename VectorType>
    class Model;

    class ModelBase;
}

template <int dim>
class Spheres
{
    // Is nonlinear solver used
    bool _isNonlinear = true; 

    // Triangulation
    std::shared_ptr<Triangulation<dim>> _triangulation;

    // Model
    std::shared_ptr<PDE::ModelWithKernels<dim, VectorType>> _model;

    // Analysis
    std::shared_ptr<CORE::AnalysisOneStep<dim, VectorType>> _analysis;

    // Create problem set
    std::shared_ptr<CORE::ProblemSet<dim, VectorType>> _problemSet;

    // Manager
    std::shared_ptr<CORE::ActionsManager<dim, VectorType>> _manager;

    // Main problem
    std::shared_ptr<CORE::Problem<dim, VectorType>> _problem;

    // Interpolator to apply initial cocnditions
    std::shared_ptr<PDE::FieldInterpolator<dim, VectorType>> _fieldInterpolator;

    double _timeStart;
    double _timeEnd;
    double _timeStepInit;
    double _timeStepMin;
    double _timeStepMax;

public:
    Spheres(dealii::ParameterHandler &param);

    void run();

private:
    std::shared_ptr<Triangulation<dim>> createMesh();
    std::shared_ptr<Triangulation<dim>> loadMesh(const std::string& filename);

    /*
    void setBoundaryCondtions(std::shared_ptr<CORE::Model<dim, VectorType>> model, Triangulation<dim>& triangulation, 
        const dealii::Point<dim>& O1, const dealii::Point<dim>& O2, const std::string& mode);
    */
};
