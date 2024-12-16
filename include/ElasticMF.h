//
// Created by Florian on 02.12.2024.
//

#ifndef LAPLACIAN_ELASTICMF_H
#define LAPLACIAN_ELASTICMF_H
// stolen from https://github.com/mwichro/SpecialTopicsFEM/blob/main/include/laplacian_mf.h
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <fstream>

#include "ElasticOperator.h"

using namespace dealii;

template <int dim, int degree, int n_component>
class ElasticMatrixFree
{
public:
  using Number     = double;
  using VectorType = dealii::LinearAlgebra::distributed::Vector<Number>;

  // OperatorType now includes n_components (dim for displacement in elasticity)
  using OperatorType = ElasticOperator::Operator<dim, degree, n_component>;

  ElasticMatrixFree(const Triangulation<dim> &tria);

  const auto &
  get_operator() const
  {
    return system_matrix;
  }

  void
  vmult_inverse(VectorType &dst, const VectorType &src) const;

  void
  initialize();

  DoFHandler<dim> dof_handler;

private:
  void
  setup_system();
  void
  solve();

  void
  output_results(const unsigned int cycle) const;

  const Triangulation<dim> &tria;

  // Change to a vector-valued finite element system
  FESystem<dim> fe;

  AffineConstraints<Number> constraints;

  std::shared_ptr<MatrixFree<dim, Number>> matrix_free_storage;
  OperatorType                             system_matrix;

  VectorType solution;
  VectorType system_rhs;

};

template <int dim, int degree, int n_component>
ElasticMatrixFree<dim, degree, n_component>::ElasticMatrixFree(
  const Triangulation<dim> &tria)
  : tria(tria)
  , fe(FE_Q<dim>(degree), dim) // Vector-valued FE system: degree per dimension
  , dof_handler(tria)
{}

template <int dim, int degree, int n_component>
void
ElasticMatrixFree<dim, degree, n_component>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  constraints.reinit();
  constraints.close();

  typename MatrixFree<dim, Number>::AdditionalData additional_data;
  additional_data.tasks_parallel_scheme =
    MatrixFree<dim, Number>::AdditionalData::none;
  additional_data.mapping_update_flags = update_gradients | update_JxW_values;

  matrix_free_storage = std::make_shared<MatrixFree<dim, Number>>();

  MappingQ1<dim> mapping;
  matrix_free_storage->reinit(mapping,
                              dof_handler,
                              constraints,
                              QGauss<1>(fe.degree + 1),
                              additional_data);

  matrix_free_storage->initialize_dof_vector(solution);
  matrix_free_storage->initialize_dof_vector(system_rhs);

  system_matrix.initialize(matrix_free_storage);
}




template <int dim, int degree, int n_component>
void
ElasticMatrixFree<dim, degree, n_component>::initialize()
{

  setup_system();
}

template <int dim, int degree, int n_component>
void
ElasticMatrixFree<dim, degree, n_component>::output_results(const unsigned int cycle) const
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  std::vector<std::string> solution_names;
  for (unsigned int d = 0; d < dim; ++d)
    solution_names.emplace_back("displacement_" + std::to_string(d));

  data_out.add_data_vector(solution, solution_names);
  data_out.build_patches();

  std::ofstream output("solution-" + std::to_string(cycle) + ".vtu");
  data_out.write_vtu(output);
}
#endif // LAPLACIAN_LAPLACIANMF_H
