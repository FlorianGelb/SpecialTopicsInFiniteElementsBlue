//
// Created by Florian on 25.11.2024.
//

#include <deal.II/base/parameter_handler.h>

#ifndef LAPLACIAN_FESYSTEM_H
#define LAPLACIAN_FESYSTEM_H
using namespace dealii;
struct FESystem
{
  unsigned int poly_degree;
  unsigned int quad_order;

  static void declare_parameters(ParameterHandler &prm);

  void parse_parameters(ParameterHandler &prm);
};
#endif // LAPLACIAN_FESYSTEM_H
