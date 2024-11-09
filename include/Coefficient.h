//
// Created by Florian on 09.11.2024.
//

#ifndef LAPLACIAN_COEFFICIENT_H
#define LAPLACIAN_COEFFICIENT_H
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/vectorization.h>  // for VectorizedArray

using namespace dealii;

template <int dim>
class Coefficient : public Function<dim>
{
public:
  // Override for `value` with `double` return type
  virtual double value(const Point<dim> &p, const unsigned int component = 0) const override;

  // Templated `value` function for different types
  template <typename number>
  number value(const Point<dim, number> &p, const unsigned int component = 0) const;
};

// Templated member function definition
template <int dim>
template <typename number>
number Coefficient<dim>::value(const Point<dim, number> &p, const unsigned int /*component*/) const
{
  return 1. / (0.05 + 2. * p.square());
}

// Override for `value` with `double`
template <int dim>
double Coefficient<dim>::value(const Point<dim> &p, const unsigned int component) const
{
  return value<double>(p, component);
}

#endif // LAPLACIAN_COEFFICIENT_H
