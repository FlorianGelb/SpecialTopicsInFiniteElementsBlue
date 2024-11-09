//
// Created by Florian on 27.10.2024.
//

#include "TestCaseK.h"
#include <cmath>

template<typename T, int dim>
T TestCaseK<T, dim>::evaluate(const dealii::Point<dim> &p)
{
  float x = p[0];
  float y = p[1];
  return  1 + x + y;
}

template class TestCaseK<double, 2>;
template class TestCaseK<double, 3>;