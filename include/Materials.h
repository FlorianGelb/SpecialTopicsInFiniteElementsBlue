//
// Created by Florian on 25.11.2024.
//
#include <deal.II/base/parameter_handler.h>
#ifndef LAPLACIAN_MATERIALS_H
#define LAPLACIAN_MATERIALS_H

struct Materials
{
  double nu;
  double mu;

  static void declare_parameters(ParameterHandler &prm);

  void parse_parameters(ParameterHandler &prm);
};



#endif // LAPLACIAN_MATERIALS_H
