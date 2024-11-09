// ---------------------------------------------------------------------
//
// Copyright (C) 2024 by Luca Heltai
//
// This file is part of the bare-dealii-app application, based on the 
// deal.II library.
//
// The bare-dealii-app application is free software; you can use it, 
// redistribute it, and/or modify it under the terms of the Apache-2.0 License
// WITH LLVM-exception as published by the Free Software Foundation; either 
// version 3.0 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md
// at the top level of the bare-dealii-app distribution. 
// 
// ---------------------------------------------------------------------

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "ConstFunction.h"
#include "LaplaceSolver.h"
#include "SinFunction.h"
#include "TestCaseK.h"
#include "TestCaseK.h"
#include "TestCaseRHS.h"
#include "MatrixFreeLaplaceOperator.h"

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/mg_tools.h>

#include <deal.II/multigrid/mg_smoother.h>

#include <deal.II/matrix_free/matrix_free.h>


int
main()
{
  Triangulation<DEAL_DIMENSION> triangulation(Triangulation<2>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);
  LaplaceOperator<DEAL_DIMENSION, 1, double> system_matrix;
  MGLevelObject<LaplaceOperator<DEAL_DIMENSION, 1, float>> mg_matrices;
  MappingQ1<DEAL_DIMENSION> mapping;
  AffineConstraints<double> constraints;
  LinearAlgebra::distributed::Vector<double> solution;
  LinearAlgebra::distributed::Vector<double> dst;
  MGConstrainedDoFs mg_constrained_dofs;
  const dealii::FE_Q<DEAL_DIMENSION> fe = dealii::FE_Q<DEAL_DIMENSION>(1);
  DoFHandler<DEAL_DIMENSION> dof_handler = DoFHandler(triangulation);




  try
    {
      // TestCaseK <double, DEAL_DIMENSION> k;
      // TestCaseRHS <double, DEAL_DIMENSION> rhs;
      // LaplaceSolver<DEAL_DIMENSION> laplace_solver(triangulation, &rhs, &k);
      // laplace_solver.run("sin");




      system_matrix.clear();
      mg_matrices.clear_elements();

      dof_handler.distribute_dofs(fe);
      dof_handler.distribute_mg_dofs();
      mg_constrained_dofs.initialize(dof_handler);

      cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;


  {
    {
      typename MatrixFree<DEAL_DIMENSION, double>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
        MatrixFree<DEAL_DIMENSION, double>::AdditionalData::none;
      additional_data.mapping_update_flags =
        (update_gradients | update_JxW_values | update_quadrature_points);
      std::shared_ptr<MatrixFree<DEAL_DIMENSION, double>> system_mf_storage(
        new MatrixFree<DEAL_DIMENSION, double>());
      system_mf_storage->reinit(mapping,
                                dof_handler,
                                constraints,
                                QGauss<1>(fe.degree + 1),
                                additional_data);
      system_matrix.initialize(system_mf_storage);
    }

    system_matrix.evaluate_coefficient(Coefficient<DEAL_DIMENSION>());

    system_matrix.initialize_dof_vector(solution);
    system_matrix.initialize_dof_vector(dst);

    for (typename LinearAlgebra::distributed::Vector<double>::size_type i = 0; i < solution.locally_owned_size(); ++i)
      {
        solution[i] = 1.0;
      }

    // Compress the vector to ensure synchronization of ghost values across processors
    solution.compress(VectorOperation::insert);



  }

  {
    const unsigned int nlevels = triangulation.n_global_levels();
    mg_matrices.resize(0, nlevels - 1);

    const std::set<types::boundary_id> dirichlet_boundary_ids = {0};
    AffineConstraints<double> level_constraints;

    for (unsigned int level = 0; level < nlevels; ++level)
      {

        typename MatrixFree<DEAL_DIMENSION, float>::AdditionalData additional_data;
        additional_data.tasks_parallel_scheme =
          MatrixFree<DEAL_DIMENSION, float>::AdditionalData::none;
        additional_data.mapping_update_flags =
          (update_gradients | update_JxW_values | update_quadrature_points);
        additional_data.mg_level = level;
        std::shared_ptr<MatrixFree<DEAL_DIMENSION, float>> mg_mf_storage_level =
          std::make_shared<MatrixFree<DEAL_DIMENSION, float>>();
        mg_mf_storage_level->reinit(mapping,
                                    dof_handler,
                                    level_constraints,
                                    QGauss<1>(fe.degree + 1),
                                    additional_data);

        mg_matrices[level].initialize(mg_mf_storage_level,
                                      mg_constrained_dofs,
                                      level);
        mg_matrices[level].evaluate_coefficient(Coefficient<DEAL_DIMENSION>());
      }
  }
  std::cout << "Before vmult, solution norm squared = " << solution.norm_sqr() << std::endl;
  Coefficient<DEAL_DIMENSION> coefficient; // Create an object of Coefficient
  // Debug the coefficient values (print a few points)
  for (const auto &cell : dof_handler.active_cell_iterators()) {
      // Loop over the vertices of each cell
      for (unsigned int v = 0; v < dealii::GeometryInfo<DEAL_DIMENSION>::vertices_per_cell; ++v) {
          dealii::Point<DEAL_DIMENSION> p = cell->vertex(v);  // Get vertex coordinates for each vertex in the cell
          std::cout << "Coefficient at point " << p << " = "
                    << coefficient.value(p) << std::endl;
        }
    }


  // Perform the matrix-free operation
  system_matrix.vmult(dst, solution);


  // Ensure synchronization of the solution vector (including ghost values)
  dst.compress(VectorOperation::insert);

  // Debug: Check the norm after vmult
  std::cout << "After vmult, solution norm squared = " << dst.norm_sqr() << std::endl;

  // Print a few values of the solution
  for (unsigned int i = 0; i < dst.size(); ++i) {
      std::cout << "solution[" << i << "] = " << dst[i] << std::endl;
    }



}

  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
