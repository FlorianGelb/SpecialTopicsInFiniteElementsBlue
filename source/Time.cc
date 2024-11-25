//
// Created by Florian on 25.11.2024.
//

#include "Time.h"
#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/block_vector.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition_selector.h>
#include <deal.II/lac/solver_cg.h>


using namespace dealii;
void Time::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Time");
  {
    prm.declare_entry("End time", "1", Patterns::Double(), "End time");

    prm.declare_entry("Time step size",
                      "0.1",
                      Patterns::Double(),
                      "Time step size");
  }
  prm.leave_subsection();
}

void Time::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Time");
  {
    end_time = prm.get_double("End time");
    delta_t  = prm.get_double("Time step size");
  }
  prm.leave_subsection();
}