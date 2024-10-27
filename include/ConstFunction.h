//
// Created by Florian on 27.10.2024.
//

#ifndef LAPLACIAN_CONSTFUNCTION_H
#define LAPLACIAN_CONSTFUNCTION_H
#include "BaseFunction.h"

template<typename T, int dim>
class ConstFunction : public BaseFunction<T, dim>
{
public:
  T evaluate(const dealii::Point<dim> &p) override;
  void setConstant(double constant);

private:
  double c;
};



#endif // LAPLACIAN_CONSTFUNCTION_H
