#include "Coefficient.h"

// Explicit instantiations for Coefficient<2> with different types for the `value` function

template class Coefficient<2>;  // Instantiates Coefficient<2> with default double value function

// Explicit instantiation of `value` for VectorizedArray<double, 1ul> and VectorizedArray<float, 1ul>
template dealii::VectorizedArray<double, 1> Coefficient<2>::value<dealii::VectorizedArray<double, 1>>(const Point<2, dealii::VectorizedArray<double, 1>> &, const unsigned int) const;

template dealii::VectorizedArray<float, 1> Coefficient<2>::value<dealii::VectorizedArray<float, 1>>(const Point<2, dealii::VectorizedArray<float, 1>> &, const unsigned int) const;