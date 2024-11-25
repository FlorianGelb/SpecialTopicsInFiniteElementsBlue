//
// Created by Florian on 25.11.2024.
//

#ifndef LAPLACIAN_POINTHISTORY_H
#define LAPLACIAN_POINTHISTORY_H


#include <deal.II/physics/elasticity/kinematics.h>
#include <deal.II/physics/elasticity/standard_tensors.h>
#include "Material_Compressible_Neo_Hook_Three_Field.h"
#include "AllParameters.h"

using namespace dealii;
template <int dim>
class PointHistory
{
public:
  PointHistory()
    : F_inv(Physics::Elasticity::StandardTensors<dim>::I)
    , tau(SymmetricTensor<2, dim>())
    , d2Psi_vol_dJ2(0.0)
    , dPsi_vol_dJ(0.0)
    , Jc(SymmetricTensor<4, dim>())
  {}

  virtual ~PointHistory() = default;

  void setup_lqp(const AllParameters &parameters)
  {
    material =
      std::make_shared<Material_Compressible_Neo_Hook_Three_Field<dim>>(
        parameters.mu, parameters.nu);
    update_values(Tensor<2, dim>(), 0.0, 1.0);
  }

  void update_values(const Tensor<2, dim> &Grad_u_n,
                const double          p_tilde,
                const double          J_tilde)
  {
    const Tensor<2, dim> F = Physics::Elasticity::Kinematics::F(Grad_u_n);
    material->update_material_data(F, p_tilde, J_tilde);

    F_inv         = invert(F);
    tau           = material->get_tau();
    Jc            = material->get_Jc();
    dPsi_vol_dJ   = material->get_dPsi_vol_dJ();
    d2Psi_vol_dJ2 = material->get_d2Psi_vol_dJ2();
  }

  double get_J_tilde() const
  {
    return material->get_J_tilde();
  }

  double get_det_F() const
  {
    return material->get_det_F();
  }

  const Tensor<2, dim> &get_F_inv() const
  {
    return F_inv;
  }

  double get_p_tilde() const
  {
    return material->get_p_tilde();
  }

  const SymmetricTensor<2, dim> &get_tau() const
  {
    return tau;
  }

  double get_dPsi_vol_dJ() const
  {
    return dPsi_vol_dJ;
  }

  double get_d2Psi_vol_dJ2() const
  {
    return d2Psi_vol_dJ2;
  }

  const SymmetricTensor<4, dim> &get_Jc() const
  {
    return Jc;
  }

private:
  std::shared_ptr<Material_Compressible_Neo_Hook_Three_Field<dim>> material;

  Tensor<2, dim> F_inv;

  SymmetricTensor<2, dim> tau;
  double                  d2Psi_vol_dJ2;
  double                  dPsi_vol_dJ;

  SymmetricTensor<4, dim> Jc;
};



#endif // LAPLACIAN_POINTHISTORY_H
