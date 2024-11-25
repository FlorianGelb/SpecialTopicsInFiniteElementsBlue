//
// Created by Florian on 25.11.2024.
//

#include "Geometry.h"
void Geometry::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Geometry");
  {
    prm.declare_entry("Global refinement",
                      "2",
                      Patterns::Integer(0),
                      "Global refinement level");

    prm.declare_entry("Grid scale",
                      "1e-3",
                      Patterns::Double(0.0),
                      "Global grid scaling factor");

    prm.declare_entry("Pressure ratio p/p0",
                      "100",
                      Patterns::Selection("20|40|60|80|100"),
                      "Ratio of applied pressure to reference pressure");
  }
  prm.leave_subsection();
}

void Geometry::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Geometry");
  {
    global_refinement = prm.get_integer("Global refinement");
    scale             = prm.get_double("Grid scale");
    p_p0              = prm.get_double("Pressure ratio p/p0");
  }
  prm.leave_subsection();
}