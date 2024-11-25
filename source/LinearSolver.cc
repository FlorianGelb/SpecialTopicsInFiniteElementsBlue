//
// Created by Florian on 25.11.2024.
//

#include "LinearSolver.h"

namespace dealii
{

  void LinearSolver::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Linear solver");
    {
      prm.declare_entry("Solver type",
                        "CG",
                        Patterns::Selection("CG|Direct"),
                        "Type of solver used to solve the linear system");

      prm.declare_entry("Residual",
                        "1e-6",
                        Patterns::Double(0.0),
                        "Linear solver residual (scaled by residual norm)");

      prm.declare_entry(
        "Max iteration multiplier",
        "1",
        Patterns::Double(0.0),
        "Linear solver iterations (multiples of the system matrix size)");

      prm.declare_entry("Use static condensation",
                        "true",
                        Patterns::Bool(),
                        "Solve the full block system or a reduced problem");

      prm.declare_entry("Preconditioner type",
                        "ssor",
                        Patterns::Selection("jacobi|ssor"),
                        "Type of preconditioner");

      prm.declare_entry("Preconditioner relaxation",
                        "0.65",
                        Patterns::Double(0.0),
                        "Preconditioner relaxation value");
    }
    prm.leave_subsection();
  }

  void LinearSolver::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Linear solver");
    {
      type_lin                  = prm.get("Solver type");
      tol_lin                   = prm.get_double("Residual");
      max_iterations_lin        = prm.get_double("Max iteration multiplier");
      use_static_condensation   = prm.get_bool("Use static condensation");
      preconditioner_type       = prm.get("Preconditioner type");
      preconditioner_relaxation = prm.get_double("Preconditioner relaxation");
    }
    prm.leave_subsection();
  }


} // namespace dealii