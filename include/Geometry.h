//
// Created by Florian on 25.11.2024.
//
#include <deal.II/base/parameter_handler.h>
#ifndef LAPLACIAN_GEOMETRY_H
#define LAPLACIAN_GEOMETRY_H


using namespace dealii;
struct Geometry
{
  unsigned int global_refinement;
  double       scale;
  double       p_p0;

  static void declare_parameters(ParameterHandler &prm);

  void parse_parameters(ParameterHandler &prm);
};



#endif // LAPLACIAN_GEOMETRY_H
