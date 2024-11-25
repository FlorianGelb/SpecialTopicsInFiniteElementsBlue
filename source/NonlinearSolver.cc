//
// Created by Florian on 25.11.2024.
//

#include "NonlinearSolver.h"
using namespace dealii;
void NonlinearSolver::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Nonlinear solver");
  {
    max_iterations_NR = prm.get_integer("Max iterations Newton-Raphson");
    tol_f             = prm.get_double("Tolerance force");
    tol_u             = prm.get_double("Tolerance displacement");
  }
  prm.leave_subsection();
}
void NonlinearSolver::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Nonlinear solver");
  {
    prm.declare_entry("Max iterations Newton-Raphson",
                      "10",
                      Patterns::Integer(0),
                      "Number of Newton-Raphson iterations allowed");

    prm.declare_entry("Tolerance force",
                      "1.0e-9",
                      Patterns::Double(0.0),
                      "Force residual tolerance");

    prm.declare_entry("Tolerance displacement",
                      "1.0e-6",
                      Patterns::Double(0.0),
                      "Displacement error tolerance");
  }
  prm.leave_subsection();
}
