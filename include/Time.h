//
// Created by Florian on 25.11.2024.
//
#ifndef LAPLACIAN_TIME_H
#define LAPLACIAN_TIME_H
#include <deal.II/base/parameter_handler.h> // Ensure this include is present
using namespace dealii;
struct Time
{
  double delta_t;
  double end_time;

  static void declare_parameters(ParameterHandler &prm);
  void parse_parameters(ParameterHandler &prm);
};



#endif // LAPLACIAN_TIME_H
