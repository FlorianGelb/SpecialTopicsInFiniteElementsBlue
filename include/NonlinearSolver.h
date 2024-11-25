//
// Created by Florian on 25.11.2024.
//

#ifndef LAPLACIAN_NONLINEARSOLVER_H
#define LAPLACIAN_NONLINEARSOLVER_H
#include <deal.II/base/parameter_handler.h>

using namespace dealii;
struct NonlinearSolver
{
  unsigned int max_iterations_NR;
  double       tol_f;
  double       tol_u;

  static void declare_parameters(ParameterHandler &prm);

  void parse_parameters(ParameterHandler &prm);
};



#endif // LAPLACIAN_NONLINEARSOLVER_H
