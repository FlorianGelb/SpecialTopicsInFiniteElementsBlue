//
// Created by Florian on 25.11.2024.
//

#ifndef LAPLACIAN_SOLID_H
#define LAPLACIAN_SOLID_H
#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/quadrature_point_data.h>
#include "TimeHandler.h"
#include "AllParameters.h"
#include "PointHistory.h"
#include "FESystem.h"

using namespace dealii;
template <int dim>
class Solid
{
public:
  Solid(const std::string &input_file);

  void run();

private:
  struct PerTaskData_ASM;
  struct ScratchData_ASM;

  struct PerTaskData_SC;
  struct ScratchData_SC;

  struct PerTaskData_UQPH;
  struct ScratchData_UQPH;

  void make_grid();

  void system_setup();

  void determine_component_extractors();

  void make_constraints(const int it_nr);

  void assemble_system();

  void assemble_system_one_cell(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    ScratchData_ASM                                      &scratch,
    PerTaskData_ASM                                      &data) const;

  void assemble_sc();

  void assemble_sc_one_cell(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    ScratchData_SC                                       &scratch,
    PerTaskData_SC                                       &data);

  void copy_local_to_global_sc(const PerTaskData_SC &data);

  void setup_qph();

  void update_qph_incremental(const BlockVector<double> &solution_delta);

  void update_qph_incremental_one_cell(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    ScratchData_UQPH                                     &scratch,
    PerTaskData_UQPH                                     &data);

  void copy_local_to_global_UQPH(const PerTaskData_UQPH & /*data*/)
  {}

  void solve_nonlinear_timestep(BlockVector<double> &solution_delta);

  std::pair<unsigned int, double>
  solve_linear_system(BlockVector<double> &newton_update);

  BlockVector<double>
  get_total_solution(const BlockVector<double> &solution_delta) const;

  void output_results() const;

  AllParameters parameters;

  double vol_reference;

  Triangulation<dim> triangulation;

  TimeHandler                time;
  mutable TimerOutput timer;

  CellDataStorage<typename Triangulation<dim>::cell_iterator,
                  PointHistory<dim>>
    quadrature_point_history;

  const unsigned int               degree;
  const dealii::FESystem<dim>              fe;
  DoFHandler<dim>                  dof_handler;
  const unsigned int               dofs_per_cell;
  const FEValuesExtractors::Vector u_fe;
  const FEValuesExtractors::Scalar p_fe;
  const FEValuesExtractors::Scalar J_fe;

  static const unsigned int n_blocks          = 3;
  static const unsigned int n_components      = dim + 2;
  static const unsigned int first_u_component = 0;
  static const unsigned int p_component       = dim;
  static const unsigned int J_component       = dim + 1;

  enum
  {
    u_dof = 0,
    p_dof = 1,
    J_dof = 2
  };

  std::vector<types::global_dof_index> dofs_per_block;
  std::vector<types::global_dof_index> element_indices_u;
  std::vector<types::global_dof_index> element_indices_p;
  std::vector<types::global_dof_index> element_indices_J;

  const QGauss<dim>     qf_cell;
  const QGauss<dim - 1> qf_face;
  const unsigned int    n_q_points;
  const unsigned int    n_q_points_f;

  AffineConstraints<double> constraints;
  BlockSparsityPattern      sparsity_pattern;
  BlockSparseMatrix<double> tangent_matrix;
  BlockVector<double>       system_rhs;
  BlockVector<double>       solution_n;

  struct Errors
  {
    Errors()
      : norm(1.0)
      , u(1.0)
      , p(1.0)
      , J(1.0)
    {}

    void reset()
    {
      norm = 1.0;
      u    = 1.0;
      p    = 1.0;
      J    = 1.0;
    }
    void normalize(const Errors &rhs)
    {
      if (rhs.norm != 0.0)
        norm /= rhs.norm;
      if (rhs.u != 0.0)
        u /= rhs.u;
      if (rhs.p != 0.0)
        p /= rhs.p;
      if (rhs.J != 0.0)
        J /= rhs.J;
    }

    double norm, u, p, J;
  };

  Errors error_residual, error_residual_0, error_residual_norm, error_update,
    error_update_0, error_update_norm;

  void get_error_residual(Errors &error_residual);

  void get_error_update(const BlockVector<double> &newton_update,
                   Errors                    &error_update);

  std::pair<double, double> get_error_dilation() const;

  double compute_vol_current() const;

  static void print_conv_header();

  void print_conv_footer();
};




#endif // LAPLACIAN_SOLID_H
