//
// Created by Florian on 27.10.2024.
//
#ifndef LAPLACIAN_TESTCASEK_H
#define LAPLACIAN_TESTCASEK_H
#include "BaseFunction.h"
#include <deal.II/base/point.h>
template <typename T, int dim>
class TestCaseK : public BaseFunction<T, dim>
{
public:
  // Override the evaluate function from the BaseFunction interface
  T evaluate(const dealii::Point<dim> &p) override;
};
#endif // LAPLACIAN_SINFUNCTION_H
