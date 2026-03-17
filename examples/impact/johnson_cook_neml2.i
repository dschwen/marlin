# Johnson-Cook rate-dependent plasticity model for NEML2
#
# Uses the inverted Johnson-Cook flow rule inside a backward-Euler update.
# A rate-independent return-mapping predictor is used as the Newton initial
# guess so the local solve starts near the plastic branch when the trial stress
# is far above yield, while the original elastic branch is preserved.
#
# Material parameters are for OFHC Copper (Cu):
# Source: Appl. Sci. 2020, 10, 2423; doi:10.3390/app10072423

[Solvers]
  [newton]
    type = NewtonWithLineSearch
    linear_solver = 'lu'
    abs_tol = 1e-8
    rel_tol = 1e-6
    max_its = 50
    verbose = true
  []
  [lu]
    type = DenseLU
  []
[]

[EquationSystems]
  [eq_sys]
    type = NonlinearSystem
    model = 'rate'
  []
[]

[Models]
  ###############################################################################
  # Compute the invariant plastic flow direction since we are doing J2 radial return
  # These quantities are computed from the trial state and remain fixed during
  # the return mapping.
  ###############################################################################
  [trial_elastic_strain]
    type = SR2LinearCombination
    to_var = 'state/Ee'
    from_var = 'forces/E old_state/Ep'
    coefficients = '1 -1'
  []
  [cauchy_stress]
    type = LinearIsotropicElasticity
    coefficient_types = 'YOUNGS_MODULUS POISSONS_RATIO'
    coefficients = '70e9 0.28'
    strain = 'state/Ee'
    stress = 'state/S'
  []
  [flow_direction]
    type = AssociativeJ2FlowDirection
    mandel_stress = 'state/S'
    flow_direction = 'forces/N'
  []
  [trial_state]
    type = ComposedModel
    models = 'trial_elastic_strain cauchy_stress flow_direction'
  []

  ###############################################################################
  # Stress update during radial return
  ###############################################################################
  [ep_rate]
    type = ScalarVariableRate
    variable = 'state/ep'
  []
  [plastic_strain_rate]
    type = AssociativePlasticFlow
    flow_direction = 'forces/N'
    flow_rate = 'state/ep_rate'
    plastic_strain_rate = 'state/Ep_rate'
  []
  [plastic_strain]
    type = SR2ForwardEulerTimeIntegration
    variable = 'state/Ep'
  []
  [plastic_update]
    type = ComposedModel
    models = 'ep_rate plastic_strain_rate plastic_strain'
  []
  [elastic_strain]
    type = SR2LinearCombination
    to_var = 'state/Ee'
    from_var = 'forces/E state/Ep'
    coefficients = '1 -1'
  []
  [stress_update]
    type = ComposedModel
    models = 'elastic_strain cauchy_stress'
  []

  ###############################################################################
  # Johnson-Cook flow rate computation
  # Uses the inverted J-C formula to compute strain rate from stress
  ###############################################################################
  [vonmises]
    type = SR2Invariant
    invariant_type = 'VONMISES'
    tensor = 'state/S'
    invariant = 'state/s'
  []
  [jc_flowrate]
    type = JohnsonCookFlowRate
    vonmises_stress = 'state/s'
    equivalent_plastic_strain = 'state/ep'
    temperature = 'forces/T'
    use_temperature = true
    flow_rate = state/ep_rate
    # OFHC Copper parameters
    A = 99.7e6
    B = 262.8e6
    n = 0.23
    C = 0.029
    m = 0.98
    reference_strain_rate = 1.0
    reference_temperature = 300
    melting_temperature = 1338
  []
  [integrate_ep]
    type = ScalarBackwardEulerTimeIntegration
    variable = 'state/ep'
  []

  ###############################################################################
  # Implicit rate model for radial return
  ###############################################################################
  [rate]
    type = ComposedModel
    models = "plastic_update stress_update vonmises jc_flowrate integrate_ep"
  []
  [radial_return]
    type = ImplicitUpdate
    equation_system = 'eq_sys'
    solver = 'newton'
    initial_guess_model = 'ep_predictor'
  []

  ###############################################################################
  # Rate-independent return-mapping predictor for Newton warm start.
  # Solves R_RI(dep) = s_trial - 3G*dep - sigma_y(ep_old + dep) = 0 via Newton
  # to get an accurate initial guess (above yield, near the actual solution).
  ###############################################################################
  [ep_predictor]
    type = JohnsonCookEPPredictor
    # inputs use default names: forces/E, old_state/Ep, old_state/ep -> state/ep
    A = 99.7e6
    B = 262.8e6
    n = 0.23
    youngs_modulus = 70e9
    poissons_ratio = 0.28
  []

  ###############################################################################
  # Put the models together (this is called by MOOSE/Bison)
  ###############################################################################
  [model]
    type = ComposedModel
    models = 'trial_state radial_return plastic_update stress_update'
    additional_outputs = 'state/Ep state/ep'
  []
[]
