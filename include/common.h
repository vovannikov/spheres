#ifndef COMMON_H
#define COMMON_H

#include <fstream>
#include <iostream>
#include <cmath>
#include <type_traits>

#include <deal.II/base/function.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/mpi.h>

#include "core/model.h"
#include "core/analysis_output.h"
#include "core/output_manager.h"
#include "core/problem.h"
#include "core/problem_set.h"
#include "core/actions_manager.h"
#include "core/util.h"

#include "pde/model_with_kernels.h"
#include "pde/postprocessor.h"

#include "tools/global.h"

#if defined(DEAL_II_WITH_TRILINOS)

#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include "core/linear_solver_trilinos.h"

//using VectorType = dealii::LinearAlgebra::distributed::Vector<double>;
using VectorType = dealii::TrilinosWrappers::MPI::Vector;
using MatrixType = dealii::TrilinosWrappers::SparseMatrix;

template<typename VectorType, typename MatrixType>
using LinearSolverType = CORE::LinearSolverTrilinos<VectorType, MatrixType>;

#else

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include "core/linear_solver_gmres.h"

using MatrixType = dealii::SparseMatrix<double>;
using VectorType = dealii::Vector<double>;

template<typename VectorType, typename MatrixType>
using LinearSolverType = CORE::LinearSolverGMRES<VectorType, MatrixType>;

#endif


#include <deal.II/distributed/shared_tria.h>
#include <deal.II/numerics/solution_transfer.h>

template<int dim>
using Triangulation = dealii::parallel::shared::Triangulation<dim>;

template<int dim>
using SolutionTransfer = dealii::SolutionTransfer<dim, VectorType>;

template<int dim>
std::shared_ptr<Triangulation<dim>> createTriaPtr()
{
    return std::make_shared<Triangulation<dim>>(TOOLS::mpi_comm);
}

#endif /* COMMON_H */