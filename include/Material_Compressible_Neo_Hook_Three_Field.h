//
// Created by Florian on 25.11.2024.
//

#ifndef LAPLACIAN_MATERIAL_COMPRESSIBLE_NEO_HOOK_THREE_FIELD_H
#define LAPLACIAN_MATERIAL_COMPRESSIBLE_NEO_HOOK_THREE_FIELD_H


#include <deal.II/physics/elasticity/kinematics.h>
#include <deal.II/physics/elasticity/standard_tensors.h>

using namespace dealii;
template <int dim>
class Material_Compressible_Neo_Hook_Three_Field
{
public:
  Material_Compressible_Neo_Hook_Three_Field(const double mu, const double nu)
    : kappa((2.0 * mu * (1.0 + nu)) / (3.0 * (1.0 - 2.0 * nu)))
    , c_1(mu / 2.0)
    , det_F(1.0)
    , p_tilde(0.0)
    , J_tilde(1.0)
    , b_bar(Physics::Elasticity::StandardTensors<dim>::I)
  {
    Assert(kappa > 0, ExcInternalError());
  }

  void update_material_data(const Tensor<2, dim> &F,
                       const double          p_tilde_in,
                       const double          J_tilde_in)
  {
    det_F                      = determinant(F);
    const Tensor<2, dim> F_bar = Physics::Elasticity::Kinematics::F_iso(F);
    b_bar                      = Physics::Elasticity::Kinematics::b(F_bar);
    p_tilde                    = p_tilde_in;
    J_tilde                    = J_tilde_in;

    Assert(det_F > 0, ExcInternalError());
  }

  SymmetricTensor<2, dim> get_tau()
  {
    return get_tau_iso() + get_tau_vol();
  }

  SymmetricTensor<4, dim> get_Jc() const
  {
    return get_Jc_vol() + get_Jc_iso();
  }

  double get_dPsi_vol_dJ() const
  {
    return (kappa / 2.0) * (J_tilde - 1.0 / J_tilde);
  }

  double get_d2Psi_vol_dJ2() const
  {
    return ((kappa / 2.0) * (1.0 + 1.0 / (J_tilde * J_tilde)));
  }

  double get_det_F() const
  {
    return det_F;
  }

  double get_p_tilde() const
  {
    return p_tilde;
  }

  double get_J_tilde() const
  {
    return J_tilde;
  }

protected:
  const double kappa;
  const double c_1;

  double                  det_F;
  double                  p_tilde;
  double                  J_tilde;
  SymmetricTensor<2, dim> b_bar;

  SymmetricTensor<2, dim> get_tau_vol() const
  {
    return p_tilde * det_F * Physics::Elasticity::StandardTensors<dim>::I;
  }

  SymmetricTensor<2, dim> get_tau_iso() const
  {
    return Physics::Elasticity::StandardTensors<dim>::dev_P * get_tau_bar();
  }

  SymmetricTensor<2, dim> get_tau_bar() const
  {
    return 2.0 * c_1 * b_bar;
  }

  SymmetricTensor<4, dim> get_Jc_vol() const
  {
    return p_tilde * det_F *
           (Physics::Elasticity::StandardTensors<dim>::IxI -
            (2.0 * Physics::Elasticity::StandardTensors<dim>::S));
  }

  SymmetricTensor<4, dim> get_Jc_iso() const
  {
    const SymmetricTensor<2, dim> tau_bar = get_tau_bar();
    const SymmetricTensor<2, dim> tau_iso = get_tau_iso();
    const SymmetricTensor<4, dim> tau_iso_x_I =
      outer_product(tau_iso, Physics::Elasticity::StandardTensors<dim>::I);
    const SymmetricTensor<4, dim> I_x_tau_iso =
      outer_product(Physics::Elasticity::StandardTensors<dim>::I, tau_iso);
    const SymmetricTensor<4, dim> c_bar = get_c_bar();

    return (2.0 / dim) * trace(tau_bar) *
             Physics::Elasticity::StandardTensors<dim>::dev_P -
           (2.0 / dim) * (tau_iso_x_I + I_x_tau_iso) +
           Physics::Elasticity::StandardTensors<dim>::dev_P * c_bar *
             Physics::Elasticity::StandardTensors<dim>::dev_P;
  }

  SymmetricTensor<4, dim> get_c_bar() const
  {
    return SymmetricTensor<4, dim>();
  }
};


#endif // LAPLACIAN_MATERIAL_COMPRESSIBLE_NEO_HOOK_THREE_FIELD_H
