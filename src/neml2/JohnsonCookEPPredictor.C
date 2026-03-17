/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/*            Marlin, a Fourier spectral solver for MOOSE             */
/*                                                                    */
/*            Copyright 2024 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifdef NEML2_ENABLED

#include "JohnsonCookEPPredictor.h"

#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/SR2.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/macaulay.h"
#include "neml2/tensors/functions/dev.h"
#include "neml2/tensors/functions/vol.h"
#include "neml2/tensors/functions/norm.h"
#include "neml2/misc/types.h"

namespace neml2
{
register_NEML2_object(JohnsonCookEPPredictor);

OptionSet
JohnsonCookEPPredictor::expected_options()
{
  OptionSet options = Model::expected_options();
  options.doc() =
      "Provides a better initial guess for the backward-Euler Newton solve in J2 return mapping "
      "with Johnson-Cook hardening. Computes the rate-independent return-mapping estimate "
      "\\f$ \\varepsilon_p^0 = \\varepsilon_p^{old} + \\langle \\sigma_{vm}^{trial} - "
      "\\sigma_y(\\varepsilon_p^{old}) \\rangle / (3G) \\f$ "
      "where the trial von Mises stress is computed internally from the total and plastic strains. "
      "Use as the 'initial_guess_model' in ImplicitUpdate.";

  // Input variable defaults match the standard MOOSE/NEML2 material driver names
  options.set_input("total_strain") = VariableName(FORCES, "E");
  options.set("total_strain").doc() = "Total strain (forces/E)";

  options.set_input("old_plastic_strain") = VariableName(OLD_STATE, "Ep");
  options.set("old_plastic_strain").doc() = "Old plastic strain (old_state/Ep)";

  options.set_input("old_equivalent_plastic_strain") = VariableName(OLD_STATE, "ep");
  options.set("old_equivalent_plastic_strain").doc() =
      "Old equivalent plastic strain (old_state/ep)";

  options.set_output("equivalent_plastic_strain") = VariableName(STATE, "ep");
  options.set("equivalent_plastic_strain").doc() =
      "Initial guess for equivalent plastic strain (state/ep). "
      "Consumed by ImplicitUpdate to warm-start Newton.";

  // Johnson-Cook hardening parameters (must match those in JohnsonCookFlowRate)
  options.set_parameter<TensorName<Scalar>>("A");
  options.set("A").doc() = "Reference yield stress (Pa)";

  options.set_parameter<TensorName<Scalar>>("B");
  options.set("B").doc() = "Hardening coefficient (Pa)";

  options.set_parameter<TensorName<Scalar>>("n");
  options.set("n").doc() = "Strain hardening exponent";

  options.set<double>("youngs_modulus");
  options.set("youngs_modulus").doc() = "Young's modulus (Pa)";

  options.set<double>("poissons_ratio");
  options.set("poissons_ratio").doc() = "Poisson's ratio";

  options.set<double>("ep_floor") = 1e-10;
  options.set("ep_floor").doc() =
      "Minimum equivalent plastic strain for power-law evaluation (avoids 0^n singularity).";

  return options;
}

JohnsonCookEPPredictor::JohnsonCookEPPredictor(const OptionSet & options)
  : Model(options),
    _E_total(declare_input_variable<SR2>("total_strain")),
    _Ep_old(declare_input_variable<SR2>("old_plastic_strain")),
    _ep_old(declare_input_variable<Scalar>("old_equivalent_plastic_strain")),
    _ep(declare_output_variable<Scalar>("equivalent_plastic_strain")),
    _A(declare_parameter<Scalar>("A", "A")),
    _B(declare_parameter<Scalar>("B", "B")),
    _n(declare_parameter<Scalar>("n", "n")),
    _eps_min(options.get<double>("ep_floor"))
{
  const double E = options.get<double>("youngs_modulus");
  const double nu = options.get<double>("poissons_ratio");
  const double K = E / (3.0 * (1.0 - 2.0 * nu)); // bulk modulus
  const double G = E / (2.0 * (1.0 + nu));       // shear modulus
  _three_K = 3.0 * K;
  _two_G = 2.0 * G;
  _three_G = 3.0 * G;
}

void
JohnsonCookEPPredictor::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  const auto opts = _ep_old.options();

  // -----------------------------------------------------------------------
  // Step 1: Compute trial elastic strain and trial Cauchy stress.
  // Uses the same decomposition as LinearIsotropicElasticity:
  //   S = 3K * vol(Ee) + 2G * dev(Ee)
  // -----------------------------------------------------------------------
  const SR2 Ee_trial(_E_total - _Ep_old);
  const auto three_K_s = Scalar::full(_three_K, opts);
  const auto two_G_s = Scalar::full(_two_G, opts);
  const SR2 S_trial = three_K_s * neml2::vol(Ee_trial) + two_G_s * neml2::dev(Ee_trial);

  // -----------------------------------------------------------------------
  // Step 2: Compute trial von Mises stress.
  // vm = sqrt(3/2) * ||dev(S)||  (same as SR2Invariant with VONMISES type)
  // -----------------------------------------------------------------------
  const auto eps_val = machine_precision(Ee_trial.scalar_type());
  const SR2 S_dev = neml2::dev(S_trial);
  const Scalar s_trial = std::sqrt(1.5) * neml2::norm(S_dev, eps_val);

  // -----------------------------------------------------------------------
  // Step 3: Approximate yield stress at old ep.
  //   sigma_y0 = A + B * max(ep_old, eps_min)^n
  // -----------------------------------------------------------------------
  const auto eps_floor = Scalar::full(_eps_min, opts);
  const auto ep_safe = _ep_old + eps_floor;
  const auto sigma_y_0 = _A + _B * pow(ep_safe, _n);

  // -----------------------------------------------------------------------
  // Step 4: Rate-independent return estimate.
  // Rate-independent return mapping with hardening.
  // Solves R(dep) = s_trial - 3G*dep - sigma_y(ep_old + dep) = 0 for dep >= 0 via Newton.
  // Using the accurate RI return (not just ep_pred = (s_trial - A)/3G) ensures the starting
  // point lands above yield, giving the outer Newton a meaningful (non-trivial) residual.
  // -----------------------------------------------------------------------
  const auto three_G_s = Scalar::full(_three_G, opts);
  const auto one = Scalar::full(1.0, opts);

  // Initial guess: zero-hardening elastic-plastic return
  auto dep_ri = macaulay(s_trial - sigma_y_0) / three_G_s;

  // Newton loop for rate-independent return (8 iterations is always sufficient)
  constexpr int max_iter_ri = 8;
  for (int k = 0; k < max_iter_ri; ++k)
  {
    const auto ep_k = _ep_old + dep_ri;
    const auto ep_safe_k = ep_k + eps_floor;
    const auto H_k = _A + _B * pow(ep_safe_k, _n);
    const auto R_k = s_trial - three_G_s * dep_ri - H_k;
    const auto dH_k = _B * _n * pow(ep_safe_k, _n - one);
    const auto dR_k = -three_G_s - dH_k; // dR/d(dep): always negative → well-conditioned
    dep_ri = macaulay(dep_ri - R_k / dR_k);
  }

  const auto ep_pred = _ep_old + dep_ri;

  if (out)
    _ep = ep_pred;

  // Derivatives: ImplicitUpdate calls this with dout=false (value only).
  // We still fill them in case the predictor is used in an AD context.
  if (dout_din)
  {
    const auto excess = macaulay(s_trial - sigma_y_0);
    // Smooth Heaviside approximation: d(macaulay(x))/dx ≈ x/(x+eps) for x>0
    const auto above_yield = excess / (excess + eps_floor);

    // d(s_trial)/d(Ee_trial): from SR2Invariant VONMISES — derivative is 3/2 * S_dev / s_trial
    const auto ds_dEe = Scalar::full(1.5, opts) * S_dev / s_trial; // SR2

    // d(s_trial)/d(E_total) = d(s_trial)/d(Ee) * d(Ee)/d(E_total) = ds_dEe (Ee = E - Ep)
    // d(s_trial)/d(Ep_old) = -ds_dEe

    // d(ep_pred)/d(E_total) = (above_yield / 3G) * ds_dEe
    if (_E_total.is_dependent())
      _ep.d(_E_total) = (above_yield / three_G_s) * ds_dEe; // SR2

    // d(ep_pred)/d(Ep_old) = -(above_yield / 3G) * ds_dEe
    if (_Ep_old.is_dependent())
      _ep.d(_Ep_old) = -(above_yield / three_G_s) * ds_dEe; // SR2

    // d(ep_pred)/d(ep_old): direct path (1) + through sigma_y hardening (-above_yield * dH/dep /
    // 3G)
    if (_ep_old.is_dependent())
    {
      const auto dH_dep = _B * _n * pow(ep_safe, _n - one);
      _ep.d(_ep_old) = one - above_yield * dH_dep / three_G_s; // Scalar
    }
  }
}
} // namespace neml2

#endif // NEML2_ENABLED
