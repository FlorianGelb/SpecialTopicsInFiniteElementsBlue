//
// Created by Florian on 27.10.2024.
//

#include "SinFunction.h"
#include <cmath>

template<typename T, int dim>
T SinFunction<T, dim>::evaluate(const dealii::Point<dim> &p)
{
  return sin(3*p[0] + 4*p[1]);
}

template class SinFunction<double, 2>;
template class SinFunction<double, 3>;