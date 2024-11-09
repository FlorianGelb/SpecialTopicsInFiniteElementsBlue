//
// Created by Florian on 09.11.2024.
//



#ifndef LAPLACIAN_MATRIXFREELAPLACEOPERATOR_H
#define LAPLACIAN_MATRIXFREELAPLACEOPERATOR_H


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include "Coefficient.h"
#include <iostream>
using namespace dealii;

template <int dim, int fe_degree, typename number>
class LaplaceOperator
  : public MatrixFreeOperators::
      Base<dim, LinearAlgebra::distributed::Vector<number>>
{
public:
  using value_type = number;

  LaplaceOperator();

  void clear() override;

  void evaluate_coefficient(const Coefficient<dim> &coefficient_function);

  virtual void compute_diagonal() override;

private:
  virtual void apply_add(
    LinearAlgebra::distributed::Vector<number>       &dst,
    const LinearAlgebra::distributed::Vector<number> &src) const override;

  void
  local_apply(const MatrixFree<dim, number>                    &data,
              LinearAlgebra::distributed::Vector<number>       &dst,
              const LinearAlgebra::distributed::Vector<number> &src,
              const std::pair<unsigned int, unsigned int> &cell_range) const;

  void local_compute_diagonal(
    const MatrixFree<dim, number>               &data,
    LinearAlgebra::distributed::Vector<number>  &dst,
    const unsigned int                          &dummy,
    const std::pair<unsigned int, unsigned int> &cell_range) const;

  Table<2, VectorizedArray<number>> coefficient;
};

#endif // LAPLACIAN_MATRIXFREELAPLACEOPERATOR_H
