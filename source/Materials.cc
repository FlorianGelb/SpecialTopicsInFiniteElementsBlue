//
// Created by Florian on 25.11.2024.
//

#include "Materials.h"
using namespace dealii;
void Materials::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Material properties");
  {
    prm.declare_entry("Poisson's ratio",
                      "0.4999",
                      Patterns::Double(-1.0, 0.5),
                      "Poisson's ratio");

    prm.declare_entry("Shear modulus",
                      "80.194e6",
                      Patterns::Double(),
                      "Shear modulus");
  }
  prm.leave_subsection();
}

void Materials::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Material properties");
  {
    nu = prm.get_double("Poisson's ratio");
    mu = prm.get_double("Shear modulus");
  }
  prm.leave_subsection();
}