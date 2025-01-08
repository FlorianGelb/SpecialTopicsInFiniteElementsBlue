//
// Created by Florian on 02.12.2024.
//

#ifndef LAPLACIAN_ELASTICOPERATOR_H
#define LAPLACIAN_ELASTICOPERATOR_H
// stolen from https://github.com/mwichro/SpecialTopicsFEM/blob/main/include/laplace_operator.h
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include "Coefficient.h"
  using namespace dealii;

  template <int dim, int degree, int n_components, typename number>
  class ElasticOperator : public MatrixFreeOperators::Base<dim, LinearAlgebra::distributed::Vector<number>>
  {
  public:
    using Number = double;
    using VectorType = dealii::LinearAlgebra::distributed::Vector<Number>;


    //using VectorType                     = dealii::Vector<Number>;
    const static unsigned int n_q_points = degree + 1;

    ElasticOperator();

    void evaluate_coefficient(
      const Coefficient<dim> &coefficient_function);

    void
    clear() override;

    void
    compute_diagonal() override;

  protected:
    void apply_add(VectorType &dst, const VectorType &src) const override;

  private:
    void
    local_apply(const dealii::MatrixFree<dim, Number> &      data,
                VectorType &                                 dst,
                const VectorType &                           src,
                const std::pair<unsigned int, unsigned int> &cell_range) const;



    Table<2, VectorizedArray<double>> coefficient;
  };



  template <int dim, int degree, int n_components, typename number>
  void
  ElasticOperator<dim, degree, n_components, number>::apply_add(
    ElasticOperator::VectorType       &dst,
    const ElasticOperator::VectorType &src) const
  {    // Ensure 'this->data' is valid and invoke 'local_apply' in a cell-wise loop.
    this->data->cell_loop(&ElasticOperator::local_apply, this, dst, src);}

  template <int dim, int degree, int n_components, typename number>
  void ElasticOperator<dim, degree, n_components, number>::evaluate_coefficient(
    const Coefficient<dim> &coefficient_function)
  {
    const unsigned int n_cells = this->data->n_cell_batches();
    FEEvaluation<dim, degree, degree + 1, 1> phi(*this->data);

    coefficient.reinit(n_cells, phi.n_q_points);
    for (unsigned int cell = 0; cell < n_cells; ++cell)
      {
        phi.reinit(cell);
        for (const unsigned int q : phi.quadrature_point_indices())
          coefficient(cell, q) =
            coefficient_function.value(phi.quadrature_point(q));
      }
  }

  template <int dim, int degree, int n_components, typename number>
  ElasticOperator<dim, degree, n_components, number>::ElasticOperator()
    : dealii::MatrixFreeOperators::Base<dim>()
  {}

  template <int dim, int degree, int n_components, typename number>
  void
  ElasticOperator<dim, degree, n_components, number>::clear()
  {
    dealii::MatrixFreeOperators::Base<dim>::clear();
  }

  template <int dim, int degree, int n_components, typename number>
  void
  ElasticOperator<dim, degree, n_components, number>::compute_diagonal()
  {
    AssertThrow(false, ExcMessage("Not implemented"));
  }

  template <int dim, int degree, int n_components, typename number>
  void
  ElasticOperator<dim, degree, n_components, number>::local_apply(
    const dealii::MatrixFree<dim, Number> &      data,
    VectorType &                                 dst,
    const VectorType &                           src,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    dealii::FEEvaluation<dim, degree, n_q_points, n_components, Number> fe_eval(data);
    double mu = 1;
    double lambda = 1;

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        fe_eval.reinit(cell);
        fe_eval.read_dof_values(src);
        fe_eval.evaluate(EvaluationFlags::gradients);
        for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
          {
            // Symmetric gradient at quadrature point
            const auto sym_grad = fe_eval.get_symmetric_gradient(q);

            // Divergence at quadrature point
            const auto div = fe_eval.get_divergence(q);

            // Material tensor contribution
            auto stress = 2.0 * mu * sym_grad +
                          dealii::make_vectorized_array(lambda) * div * dealii::unit_symmetric_tensor<dim>();

            fe_eval.submit_symmetric_gradient(stress, q);
          }

        fe_eval.integrate(EvaluationFlags::gradients);
        fe_eval.distribute_local_to_global(dst);
      }
  }

#endif // LAPLACIAN_ELASTICOPERATOR_H
