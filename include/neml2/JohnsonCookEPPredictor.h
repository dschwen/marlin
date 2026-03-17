/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/*            Marlin, a Fourier spectral solver for MOOSE             */
/*                                                                    */
/*            Copyright 2024 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#ifdef NEML2_ENABLED

#include "neml2/models/Model.h"

namespace neml2
{
class Scalar;
class SR2;

/**
 * @brief Provides a better initial guess for the Newton solve in J2 return mapping
 *        with Johnson-Cook hardening.
 *
 * When the trial von Mises stress is far above yield, the default predictor
 * (state/ep = old_state/ep) leaves Newton with a huge residual and an ill-conditioned
 * Jacobian. This model computes the rate-independent return-mapping estimate:
 *
 *   ep_pred = ep_old + max(0, s_trial - sigma_y(ep_old)) / (3G)
 *
 * where s_trial is computed internally from the raw strain inputs, and sigma_y = A + B*ep^n.
 *
 * Use as the `initial_guess_model` in `ImplicitUpdate`. The model is evaluated once before
 * the Newton loop starts; its `state/ep` output replaces the default `old_state/ep` predictor.
 * Because the output is NOT declared in the outer ComposedModel (it is only used by
 * ImplicitUpdate internally), there is no "multiple providers" conflict.
 */
class JohnsonCookEPPredictor : public Model
{
public:
  static OptionSet expected_options();

  JohnsonCookEPPredictor(const OptionSet & options);

protected:
  void set_value(bool out, bool dout_din, bool d2out_din2) override;

  /// Total strain (forces/E)
  const Variable<SR2> & _E_total;

  /// Old plastic strain (old_state/Ep)
  const Variable<SR2> & _Ep_old;

  /// Old equivalent plastic strain (old_state/ep)
  const Variable<Scalar> & _ep_old;

  /// Equivalent plastic strain output (state/ep) — the initial guess for Newton
  Variable<Scalar> & _ep;

  /// Reference yield stress (A)
  const Scalar & _A;

  /// Hardening coefficient (B)
  const Scalar & _B;

  /// Strain hardening exponent (n)
  const Scalar & _n;

  /// 3 * bulk modulus (3K = E / (1-2*nu))
  double _three_K;

  /// 2 * shear modulus (2G = E / (1+nu))
  double _two_G;

  /// 3 * shear modulus (3G = 3E / (2*(1+nu)))
  double _three_G;

  /// Floor for ep in the power-law derivative, avoids 0^n
  double _eps_min;
};
} // namespace neml2

#endif // NEML2_ENABLED
