#ifndef LINEAR_SOLVER_SP_H
#define LINEAR_SOLVER_SP_H

#include "core/linear_solver.h"

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

#include "tools/global.h"

template <typename VectorType, typename MatrixType>
class LinearSolverSp : public CORE::LinearSolver<VectorType, MatrixType>
{
private:
    unsigned int _maxIter;
    double _absTol;
    double _relTol;

    // create preconditioner (AMG form Trilinos ML)
    dealii::TrilinosWrappers::PreconditionAMG _preconditioner;

    // No preconditioner
    //dealii::TrilinosWrappers::PreconditionIdentity _preconditioner;

    // ILU preconditioner
    //dealii::TrilinosWrappers::PreconditionILU _preconditioner;

    // ILUT preconditioner
    //dealii::TrilinosWrappers::PreconditionILUT _preconditioner;

public:
    LinearSolverSp(unsigned int maxIter, double absTol, double relTol) 
        : _maxIter(maxIter)
        , _absTol(absTol)
        , _relTol(relTol)
    {
    }

    virtual bool solve(CORE::LinearSOE<VectorType, MatrixType>& linearSOE)
    {
        dealii::TimerOutput::Scope timer_section(TOOLS::timer::instance(), "CORE::LinearSolverSp::solve()");
        bool has_solver_succeeded;

        // for ILU
        /*
        dealii::TrilinosWrappers::PreconditionILU::AdditionalData ad;
        ad.ilu_fill = 0;
        ad.ilu_atol = 0;
        ad.ilu_rtol = 1;
        ad.overlap = 0;
        _preconditioner.initialize(linearSOE.getMatrix(), ad);
        */
        
        // Default init call
        //_preconditioner.initialize(linearSOE.getMatrix());

        // for ILUT
        //_preconditioner.initialize(linearSOE.getMatrix());

        // for AMG
        
        dealii::TrilinosWrappers::PreconditionAMG::AdditionalData ad;
        ad.smoother_sweeps = 3;
        ad.n_cycles        = 4;
        ad.smoother_type   = "Chebyshev";
        ad.coarse_type     = "Amesos-KLU";
        _preconditioner.initialize(linearSOE.getMatrix(), ad);
        
        
        // use the CG solver from deal.II
        dealii::SolverControl solver_control(_maxIter, _absTol, _relTol);

        dealii::TrilinosWrappers::SolverCG solver(solver_control);
        //dealii::TrilinosWrappers::SolverGMRES solver(solver_control);

        try {
            solver.solve(
                linearSOE.getMatrix(), 
                linearSOE.getLHS(),
                linearSOE.getRHS(), 
                _preconditioner);
            TOOLS::cout::instance() << "      Linear solver: " << solver_control.last_step() << " iterations." << std::endl;
            has_solver_succeeded = true;
        }
        catch (const std::exception& e) {
            TOOLS::cout::instance() << e.what() << std::endl;
            has_solver_succeeded = false;
        }

        return has_solver_succeeded;
    }
};

#endif /* LINEAR_SOLVER_SP_H */
