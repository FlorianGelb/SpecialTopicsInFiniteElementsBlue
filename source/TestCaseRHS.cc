//
// Created by Florian on 27.10.2024.
//

#include "TestCaseRHS.h"
#include <cmath>

template<typename T, int dim>
T TestCaseRHS<T, dim>::evaluate(const dealii::Point<dim> &p)
{
  float x = p[0];
  float y = p[1];
  return 2 * pow(x, 3) + pow(x, 2) * pow(y, 2) + 3 * pow(x, 2) * y
         - pow(x, 2) +  3 * pow(y, 2) * x
         - 7*x*y - x + 2 * pow(y, 3) - pow(y, 2) - y;
}

template class TestCaseRHS<double, 2>;
template class TestCaseRHS<double, 3>;