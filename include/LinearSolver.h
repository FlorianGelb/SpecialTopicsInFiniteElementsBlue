//
// Created by Florian on 25.11.2024.
//

#ifndef LAPLACIAN_LINEARSOLVER_H
#define LAPLACIAN_LINEARSOLVER_H
#include <deal.II/base/parameter_handler.h>
namespace dealii
{

  struct LinearSolver
  {
    std::string type_lin;
    double      tol_lin;
    double      max_iterations_lin;
    bool        use_static_condensation;
    std::string preconditioner_type;
    double      preconditioner_relaxation;

    static void declare_parameters(ParameterHandler &prm);

    void parse_parameters(ParameterHandler &prm);
  };

} // namespace dealii

#endif // LAPLACIAN_LINEARSOLVER_H
