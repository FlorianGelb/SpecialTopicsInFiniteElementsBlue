//
//
// Created by Florian on 09.11.2024.
//


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



template <int dim, int fe_degree, typename number>
LaplaceOperator<dim, fe_degree, number>::LaplaceOperator()
  : MatrixFreeOperators::Base<dim,
                              LinearAlgebra::distributed::Vector<number>>()
{}



template <int dim, int fe_degree, typename number>
void LaplaceOperator<dim, fe_degree, number>::clear()
{
  coefficient.reinit(0, 0);
  MatrixFreeOperators::Base<dim, LinearAlgebra::distributed::Vector<number>>::
    clear();
}




template <int dim, int fe_degree, typename number>
void LaplaceOperator<dim, fe_degree, number>::evaluate_coefficient(
  const Coefficient<dim> &coefficient_function)
{
  const unsigned int n_cells = this->data->n_cell_batches();
  FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(*this->data);

  coefficient.reinit(n_cells, phi.n_q_points);
  for (unsigned int cell = 0; cell < n_cells; ++cell)
    {
      phi.reinit(cell);
      for (const unsigned int q : phi.quadrature_point_indices())
        coefficient(cell, q) =
          coefficient_function.value(phi.quadrature_point(q));
    }
}




template <int dim, int fe_degree, typename number>
void LaplaceOperator<dim, fe_degree, number>::local_apply(
  const MatrixFree<dim, number>                    &data,
  LinearAlgebra::distributed::Vector<number>       &dst,
  const LinearAlgebra::distributed::Vector<number> &src,
  const std::pair<unsigned int, unsigned int>      &cell_range) const
{
  FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(data);

  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      AssertDimension(coefficient.size(0), data.n_cell_batches());
      AssertDimension(coefficient.size(1), phi.n_q_points);

      phi.reinit(cell);
      phi.read_dof_values(src);
      phi.evaluate(EvaluationFlags::gradients);
      for (const unsigned int q : phi.quadrature_point_indices())
        phi.submit_gradient(coefficient(cell, q) * phi.get_gradient(q), q);
      phi.integrate(EvaluationFlags::gradients);
      phi.distribute_local_to_global(dst);

    }
}



template <int dim, int fe_degree, typename number>
void LaplaceOperator<dim, fe_degree, number>::apply_add(
  LinearAlgebra::distributed::Vector<number>       &dst,
  const LinearAlgebra::distributed::Vector<number> &src) const
{
  this->data->cell_loop(&LaplaceOperator::local_apply, this, dst, src);
  cout << dst.norm_sqr();
}



template <int dim, int fe_degree, typename number>
void LaplaceOperator<dim, fe_degree, number>::compute_diagonal()
{
  this->inverse_diagonal_entries.reset(
    new DiagonalMatrix<LinearAlgebra::distributed::Vector<number>>());
  LinearAlgebra::distributed::Vector<number> &inverse_diagonal =
    this->inverse_diagonal_entries->get_vector();
  this->data->initialize_dof_vector(inverse_diagonal);
  unsigned int dummy = 0;
  this->data->cell_loop(&LaplaceOperator::local_compute_diagonal,
                        this,
                        inverse_diagonal,
                        dummy);

  this->set_constrained_entries_to_one(inverse_diagonal);

  for (unsigned int i = 0; i < inverse_diagonal.locally_owned_size(); ++i)
    {
      Assert(inverse_diagonal.local_element(i) > 0.,
             ExcMessage("No diagonal entry in a positive definite operator "
                        "should be zero"));
      inverse_diagonal.local_element(i) =
        1. / inverse_diagonal.local_element(i);
    }
}



template <int dim, int fe_degree, typename number>
void LaplaceOperator<dim, fe_degree, number>::local_compute_diagonal(
  const MatrixFree<dim, number>              &data,
  LinearAlgebra::distributed::Vector<number> &dst,
  const unsigned int &,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(data);

  AlignedVector<VectorizedArray<number>> diagonal(phi.dofs_per_cell);

  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      AssertDimension(coefficient.size(0), data.n_cell_batches());
      AssertDimension(coefficient.size(1), phi.n_q_points);

      phi.reinit(cell);
      for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < phi.dofs_per_cell; ++j)
            phi.submit_dof_value(VectorizedArray<number>(), j);
          phi.submit_dof_value(make_vectorized_array<number>(1.), i);

          phi.evaluate(EvaluationFlags::gradients);
          for (const unsigned int q : phi.quadrature_point_indices())
            phi.submit_gradient(coefficient(cell, q) * phi.get_gradient(q),
                                q);
          phi.integrate(EvaluationFlags::gradients);
          diagonal[i] = phi.get_dof_value(i);
        }
      for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
        phi.submit_dof_value(diagonal[i], i);
      phi.distribute_local_to_global(dst);
    }
}

template class LaplaceOperator<2, 1, double>;
template class LaplaceOperator<2, 1, float>;