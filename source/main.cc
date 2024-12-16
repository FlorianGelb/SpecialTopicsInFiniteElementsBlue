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
#include "ElasticMF.h"

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

using namespace dealii;

int
main()
{


  try
    {
      // Define the dimension for the problem
      constexpr int dim = 2; // or 3 for 3D problems

      // Create a triangulation and generate a hypercube grid
      Triangulation<DEAL_DIMENSION> triangulation;
      GridGenerator::hyper_cube(triangulation);
      triangulation.refine_global(3); // Globally refine the grid 3 times

      // Create the ElasticMatrixFree object
      constexpr int degree = 1; // Polynomial degree for FE_Q elements
      ElasticMatrixFree<DEAL_DIMENSION, degree, DEAL_DIMENSION> elastic_problem(triangulation);

      // Initialize the system
      elastic_problem.initialize();
      auto elastic_MF = elastic_problem.get_operator();

      using VectorType = dealii::LinearAlgebra::distributed::Vector<double>;

      // Loop over all dimensions
      for (unsigned int i = 0; i < dim; ++i)
        {
          // Create a unit vector e
          VectorType e; // Assuming e is initialized to the correct size elsewhere
          e.reinit(DEAL_DIMENSION); // Initialize size based on the operator
          e = 0; // Set all elements to zero
          e[i] = 1; // Set the ith component to 1

          // Create a result vector to hold the output
          VectorType result;
          result.reinit(DEAL_DIMENSION);

          // Apply the operator
          elastic_MF.vmult(result, e);

          // Output the result
          std::cout << "Result of operator applied to unit vector in direction " << i << ":\n";
          result.print(std::cout);
          std::cout << std::endl;
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