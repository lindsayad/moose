################################################################################
## Molten Salt Fast Reactor - Euratom EVOL + Rosatom MARS Design              ##
## Pronghorn input file to initialize velocity fields                         ##
## This runs a slow relaxation to steady state while ramping down the fluid   ##
## viscosity.                                                                 ##
################################################################################

# Material properties
rho = 3279.  # density [kg / m^3]  (@1000K)
mu = 0.005926 # viscosity [Pa s]

# Mass flow rate tuning
friction = 11.0  # [kg / m^4]
pump_force = ${fparse 0.25*4.0e6} # [N / m^3]
porosity = 1.0

# Numerical scheme parameters
advected_interp_method = 'upwind'
velocity_interp_method = 'rc'

[GlobalParams]
  rhie_chow_user_object = 'pins_rhie_chow_interpolator'

  two_term_boundary_expansion = true
  advected_interp_method = ${advected_interp_method}
  velocity_interp_method = ${velocity_interp_method}
  u = superficial_vel_x
  v = superficial_vel_y
  w = superficial_vel_z
  pressure = pressure
  porosity = porosity

  rho = rho
  mu = mu

[]

################################################################################
# GEOMETRY
################################################################################
[Mesh]
  [restart]
    type = FileMeshGenerator
    use_for_exodus_restart = true
    file = 'run_ns_0.e'
  []
[]


[Problem]
  kernel_coverage_check = false
[]

################################################################################
# EQUATIONS: VARIABLES, KERNELS & BCS
################################################################################

[UserObjects]
  [pins_rhie_chow_interpolator]
    type = PINSFVRhieChowInterpolator
  []
[]

[Variables]
  [pressure]
    type = INSFVPressureVariable
    initial_from_file_var = pressure
  []
  [superficial_vel_x]
    type = PINSFVSuperficialVelocityVariable
    initial_from_file_var = superficial_vel_x
  []
  [superficial_vel_y]
    type = PINSFVSuperficialVelocityVariable
    initial_from_file_var = superficial_vel_y
  []
  [superficial_vel_z]
    type = PINSFVSuperficialVelocityVariable
    initial_from_file_var = superficial_vel_z
  []
  #[lambda]
  #  family = SCALAR
  #  order = FIRST
  #  #initial_from_file_var = lambda
  #[]
[]

[AuxVariables]
  [eddy_viscosity]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
[]

[AuxKernels]
  [turbulent_viscosity]
    type = INSFVMixingLengthTurbulentViscosityAux
    variable = eddy_viscosity
    mixing_length = ${fparse 0.07*0.1}
  []
[]

[FVKernels]

  [mass]
    type = PINSFVMassAdvection
    variable = pressure
    advected_interp_method = 'skewness-corrected'
    velocity_interp_method = 'rc'
    rho = ${rho}
  []
  #[mean_zero_pressure]
  #  type = FVIntegralValueConstraint
  #  variable = pressure
  #  lambda = lambda
  #[]

  [u_time]
    type = PINSFVMomentumTimeDerivative
    variable = superficial_vel_x
    rho = ${rho}
    momentum_component = 'x'
  []
  [u_advection]
    type = PINSFVMomentumAdvection
    variable = superficial_vel_x
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
    momentum_component = 'x'
  []
  [u_viscosity]
    type = PINSFVMomentumDiffusion
    variable = superficial_vel_x
    momentum_component = 'x'
  []
  [u_viscosity_rans]
    type = INSFVMixingLengthReynoldsStress
    variable = superficial_vel_x
    rho = ${rho}
    mixing_length = ${fparse 0.07*0.1 * 1}
    momentum_component = 'x'
  []
  [u_pressure]
    type = PINSFVMomentumPressure
    variable = superficial_vel_x
    momentum_component = 'x'
    pressure = pressure
  []
  [u_friction_pump]
    type = PINSFVMomentumFriction
    variable = superficial_vel_x
    momentum_component = 'x'
    Darcy_name = 'DFC'
    Forchheimer_name = 'FFC'
    block = 'pump'
  []
  [u_friction_pump_correction]
    type = PINSFVMomentumFrictionCorrection
    variable = superficial_vel_x
    momentum_component = 'x'
    porosity = porosity
    rho = ${rho}
    Darcy_name = 'DFC'
    Forchheimer_name = 'FFC'
    block = 'mixing-plate'
    consistent_scaling = 100.
  []

  [u_friction_mixing_plate]
    type = PINSFVMomentumFriction
    variable = superficial_vel_x
    momentum_component = 'x'
    Darcy_name = 'DFC_plate'
    Forchheimer_name = 'FFC_plate'
    block = 'mixing-plate'
  []
  [u_friction_mixing_plate_correction]
    type = PINSFVMomentumFrictionCorrection
    variable = superficial_vel_x
    momentum_component = 'x'
    porosity = porosity
    rho = ${rho}
    Darcy_name = 'DFC_plate'
    Forchheimer_name = 'FFC_plate'
    block = 'mixing-plate'
    consistent_scaling = 100.
  []


  [v_time]
    type = PINSFVMomentumTimeDerivative
    variable = superficial_vel_y
    rho = ${rho}
    momentum_component = 'y'
  []
  [v_advection]
    type = PINSFVMomentumAdvection
    variable = superficial_vel_y
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
    momentum_component = 'y'
  []
  [v_viscosity]
    type = PINSFVMomentumDiffusion
    variable = superficial_vel_y
    momentum_component = 'y'
  []
  [v_viscosity_rans]
    type = INSFVMixingLengthReynoldsStress
    variable = superficial_vel_y
    rho = ${rho}
    mixing_length = ${fparse 0.07*0.1 * 1}
    momentum_component = 'y'
  []
  [v_pressure]
    type = PINSFVMomentumPressure
    variable = superficial_vel_y
    momentum_component = 'y'
    pressure = pressure
  []
  [v_friction_pump]
    type = PINSFVMomentumFriction
    variable = superficial_vel_y
    momentum_component = 'y'
    Darcy_name = 'DFC'
    Forchheimer_name = 'FFC'
    block = 'pump'
  []
  [v_friction_pump_correction]
    type = PINSFVMomentumFrictionCorrection
    variable = superficial_vel_y
    momentum_component = 'y'
    porosity = porosity
    rho = ${rho}
    Darcy_name = 'DFC'
    Forchheimer_name = 'FFC'
    block = 'mixing-plate'
    consistent_scaling = 100.
  []
  [v_friction_mixing_plate]
    type = PINSFVMomentumFriction
    variable = superficial_vel_y
    momentum_component = 'y'
    Darcy_name = 'DFC_plate'
    Forchheimer_name = 'FFC_plate'
    block = 'mixing-plate'
  []
  [pump]
    type = INSFVBodyForce
    variable = superficial_vel_y
    functor = ${pump_force}
    block = 'pump'
    momentum_component = 'y'
  []

  [w_time]
    type = PINSFVMomentumTimeDerivative
    variable = superficial_vel_z
    rho = ${rho}
    momentum_component = 'z'
  []
  [w_advection]
    type = PINSFVMomentumAdvection
    variable = superficial_vel_z
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
    momentum_component = 'z'
  []
  [w_viscosity]
    type = PINSFVMomentumDiffusion
    variable = superficial_vel_z
    momentum_component = 'z'
  []
  [w_viscosity_rans]
    type = INSFVMixingLengthReynoldsStress
    variable = superficial_vel_z
    rho = ${rho}
    mixing_length = ${fparse 0.07*0.1 * 1}
    momentum_component = 'z'
  []
  [w_pressure]
    type = PINSFVMomentumPressure
    variable = superficial_vel_z
    momentum_component = 'z'
    pressure = pressure
  []
  [w_friction_pump]
    type = PINSFVMomentumFriction
    variable = superficial_vel_z
    momentum_component = 'z'
    Darcy_name = 'DFC'
    Forchheimer_name = 'FFC'
    block = 'pump'
  []
  [w_friction_mixing_plate]
    type = PINSFVMomentumFriction
    variable = superficial_vel_z
    momentum_component = 'z'
    Darcy_name = 'DFC_plate'
    Forchheimer_name = 'FFC_plate'
    block = 'mixing-plate'
  []

[]

[FVBCs]
  [no-slip-u]
    type = INSFVNoSlipWallBC
    boundary = 'wall-reactor wall-pipe wall-pump wall-reactor-full'
    variable = superficial_vel_x
    function = 0
  []
  [no-slip-v]
    type = INSFVNoSlipWallBC
    boundary = 'wall-reactor wall-pipe wall-pump wall-reactor-full'
    variable = superficial_vel_y
    function = 0
  []
  [no-slip-w]
    type = INSFVNoSlipWallBC
    boundary = 'wall-reactor wall-pipe wall-pump wall-reactor-full'
    variable = superficial_vel_z
    function = 0
  []
[]

################################################################################
# MATERIALS
################################################################################

[Functions]
  [ad_rampdown_mu_func]
    type = ADParsedFunction
    expression = mu*(0.1*exp(-3*t)+1)
    symbol_names = 'mu'
    symbol_values = ${mu}
  []
  [mu_x]
    type = ADParsedFunction
    expression = 'if(x > -0.08 & x < 0.2 , 10.*mu*(0.1*exp(-3*t)+1),mu*(0.1*exp(-3*t)+1))'
    symbol_names = 'mu'
    symbol_values = ${mu}
  []
[]

[Materials]
  [generic]
    type = ADGenericFunctorMaterial      #defines mu artificially for numerical convergence
    prop_names = 'mu rho porosity'       #it converges to the real mu eventually.
    prop_values = '${mu} ${rho} ${porosity}'
  []
  #[mu_spatial]
  #  type = ADPiecewiseByBlockFunctorMaterial
  #  prop_name = 'mu'
  #  subdomain_to_prop_value = 'pipe mu
  #                             pump mu
  #                             mixing-plate mu
  #                             reactor mu'
  #[]
  [friction_material_pump]
    type = ADGenericVectorFunctorMaterial      #defines mu artificially for numerical convergence
    prop_names = 'DFC FFC'     #it converges to the real mu eventually.
    prop_values = '${fparse 1*friction} ${fparse 1*friction} ${fparse 1*friction} 
                   ${fparse 1*friction} ${fparse 1*friction} ${fparse 1*friction}'
  []
  [friction_material_mixing_plate]
    type = ADGenericVectorFunctorMaterial      #defines mu artificially for numerical convergence
    prop_names = 'DFC_plate FFC_plate'     #it converges to the real mu eventually.
    prop_values = '${fparse 1*friction} ${fparse 1*friction} ${fparse 1*friction} 
                   ${fparse 1*friction} ${fparse 1*friction} ${fparse 1*friction}'
  []
[]

################################################################################
# EXECUTION / SOLVE
################################################################################
[Preconditioning]
  active = FSP
  [FSP]
    type = FSP
    topsplit = 'up'
      [up]
        splitting = 'us p'
        splitting_type  = schur
        petsc_options = '-ksp_monitor_true_residual'
        petsc_options_iname = '-pc_fieldsplit_schur_fact_type  -pc_fieldsplit_schur_precondition -ksp_gmres_restart -ksp_rtol -ksp_type -ksp_atol'
        petsc_options_value = 'full                            self                             300                 1e-2      fgmres    3e-5'
        vars = 'superficial_vel_x superficial_vel_y superficial_vel_z pressure'
      []
        [us]
          vars = 'superficial_vel_x superficial_vel_y superficial_vel_z'
          splitting = 'u v w'
          splitting_type = symmetric_multiplicative
          petsc_options_iname = ' -ksp_gmres_restart -ksp_rtol -ksp_type -ksp_atol -ksp_converged_reason'
          petsc_options_value = '300                 1e-1      fgmres 1e-9          ::failed'
        []
          [u]
            vars = 'superficial_vel_x'
            petsc_options_iname = '-pc_type -ksp_pc_side -ksp_type -ksp_rtol'
            petsc_options_value = 'lu       right        gmres     1e-5'
          []
          [v]
            vars = 'superficial_vel_y'
            petsc_options_iname = '-pc_type -ksp_pc_side -ksp_type -ksp_rtol'
            petsc_options_value = 'lu       right        gmres     1e-5'
          []
          [w]
            vars = 'superficial_vel_z'
            petsc_options_iname = '-pc_type -ksp_pc_side -ksp_type -ksp_rtol'
            petsc_options_value = 'lu       right        gmres     1e-5'
          []
        [p]
          vars = 'pressure'
          petsc_options = '-ksp_monitor_true_residual -pc_lsc_scale_diag'
          petsc_options_iname = '-ksp_type -ksp_gmres_restart -ksp_rtol -pc_type -ksp_pc_side -lsc_pc_type'
          petsc_options_value = 'gmres     300                0.5       lsc       right       lu'
        []
    []
  [SMP]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_shift_type'
    petsc_options_value = 'lu       NONZERO'
  []
[]

[Executioner]
  type = Transient

  # Time-stepping parameters
  start_time = 0.0
  end_time = 150

  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 10
    dt = 0.5
    timestep_limiting_postprocessor = 'dt_limit'
  []

  # Solver parameters
  solve_type = 'NEWTON'
  #petsc_options_iname = '-pc_type -pc_factor_shift_type -ksp_gmres_restart'
  #petsc_options_value = 'lu NONZERO 20'
  line_search = 'none'
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-8
  nl_max_its = 10 # fail early and try again with a shorter time step
  l_max_its = 80
  automatic_scaling = true
[]

[Debug]
    show_var_residual_norms = true
[]

################################################################################
# SIMULATION OUTPUTS
################################################################################

[Outputs]
  csv = false
  hide = 'dt_limit'
  [restart]
    type = Exodus
    execute_on = 'timestep_end final'
  []
  # Reduce base output
  print_linear_converged_reason = false
  print_linear_residuals = false
  print_nonlinear_converged_reason = false
[]

[Postprocessors]
  [max_v]
    type = ElementExtremeValue
    variable = superficial_vel_x
    value_type = max
    block = 'reactor pump pipe'
  []
  [flow_hx_bot]
    type = VolumetricFlowRate
    boundary = 'reactor_bot'
    vel_x = superficial_vel_x
    vel_y = superficial_vel_y
    advected_quantity = ${rho}
  []
  [flow_hx_top]
    type = VolumetricFlowRate
    boundary = 'reactor_top'
    vel_x = superficial_vel_x
    vel_y = superficial_vel_y
    advected_quantity = ${rho}
  []
  [pdrop]
    type = PressureDrop
    pressure = pressure
    upstream_boundary = 'reactor_bot'
    downstream_boundary = 'mixing-plate-downstream'
    boundary = 'mixing-plate-downstream reactor_bot'
  []
  [dt_limit]
    type = Receiver
    default = 1
  []
[]
