//
// Created by Florian on 02.01.2025.
//

#include <gtest/gtest.h>
#include "ElasticMF.h"
#include "Elasticity.h"
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <math.h>

using namespace dealii;


TEST(ElasticOperator, VmultTest){

  Triangulation<3> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2); // Globally refine the grid 3 times

  // Create the ElasticMatrixFree object
  constexpr int degree = 1; // Polynomial degree for FE_Q elements
  ElasticMatrixFree<3, degree, 3> elastic_problem(triangulation);
  Elasticity<3> elastic_problem_2 = Elasticity(triangulation, degree);

  // Initialize the system
  elastic_problem.initialize();
  elastic_problem_2.intinlize();
  auto elastic_MF = elastic_problem.get_operator();

  auto elastic = &elastic_problem_2.system_matrix;

  using VectorType = dealii::LinearAlgebra::distributed::Vector<double>;

  // Loop over all dimensions
  float total_difference = 0;
  for (unsigned int i = 0; i < elastic_problem.dof_handler.n_dofs(); ++i)
    {

      VectorType e;
      e.reinit(elastic_problem.dof_handler.n_dofs());  // Ensure the vector has the correct size for DoF
      e = 0;  // Set all elements to zero
      e(i) = 1;  // Set the i-th component to 1 (unit vector)

      // Create a result vector to hold the output
      VectorType result;
      VectorType result_2;
      result.reinit(e.size());  // Ensure result has the same size as e
      result_2.reinit(e.size());

      // Output the initial unit vector (for debugging purposes)


      // Apply the operator
      elastic_MF.vmult(result, e);
      elastic->vmult(result_2, e);

      VectorType difference;
      difference.reinit(result.size()); // Ensure the vector is the correct size
      difference = result_2;           // Copy result_2 into the difference
      difference.add(-1.0, result);    // Compute result_2 - result

      // Output the result

      total_difference += difference.l1_norm();


    }
  EXPECT_LE(total_difference / elastic_problem.dof_handler.n_dofs(), 1E-14);

}