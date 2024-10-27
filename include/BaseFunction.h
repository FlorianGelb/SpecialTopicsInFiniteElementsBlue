//
// Created by Florian on 27.10.2024.
//

#ifndef LAPLACIAN_BASEFUNCTION_H
#define LAPLACIAN_BASEFUNCTION_H

#include <vector>
#include <deal.II/base/point.h>
template <typename T, int dim> // Enable usage of vectors in C^n or R^n
// Abstract Class / Interface to ensure evaluated(x) can be called on
// user defined classes
class BaseFunction
{
public:
  // Some C++ magic i dont know
  virtual ~BaseFunction() = default;
  // Vector is maybe not the right data type but im no C++ expert
  // and Google says arrays cannot be returned
  // (at least the first couple of results)
  virtual T evaluate(const dealii::Point<dim> &p) = 0;
};


#endif // LAPLACIAN_BASEFUNCTION_H
