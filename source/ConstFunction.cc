//
// Created by Florian on 27.10.2024.
//

#include "ConstFunction.h"

template<typename T, int dim>
void ConstFunction<T, dim>::setConstant(double constant)
{
  c = constant;
}


template<typename T, int dim>
T ConstFunction<T, dim>::evaluate(const dealii::Point<dim> &p)
{
  return c;
}

template class ConstFunction<double, 2>;
template class ConstFunction<double, 3>;