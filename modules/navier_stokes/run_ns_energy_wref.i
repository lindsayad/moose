################################################################################
## Molten Salt Fast Reactor - Euratom EVOL + Rosatom MARS Design              ##
## Pronghorn input file to initialize velocity fields                         ##
## This runs a slow relaxation to steady state while ramping down the fluid   ##
## viscosity.                                                                 ##
################################################################################

# Material properties fuel
mu = 0.005926 # viscosity [Pa s]
rho_solid = 3279.7
rho_liquid = 3279.7
k_solid = 0.38
k_liquid = 0.38
cp_solid = 640.
cp_liquid = 640.
L = 4e5
T_liquidus = 0.1 #812.
T_solidus = 0. #785.

# Material properties reflector
k_ref = 30.
cp_ref = 880.
rho_ref = 3580.

power = 10e3

# Mass flow rate tuning
friction = 11.0 # [kg / m^4]
pump_force = '${fparse 0.25*4.0e6}' # [N / m^3]
porosity = 1.0
Reactor_Area = ${fparse 3.14159*0.2575*0.2575}
Pr = ${fparse mu*cp_liquid/k_solid}

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

  rho = ${rho_liquid}
  mu = mu
  cp = 'cp_mixture'
[]

################################################################################
# GEOMETRY
################################################################################
[Mesh]
  [fmg]
    type = FileMeshGenerator
    use_for_exodus_restart = true
    file = run_ns_energy_wref_advected_final.e
  []
[]


################################################################################
# EQUATIONS: VARIABLES, KERNELS & BCS
################################################################################

[UserObjects]
  [pins_rhie_chow_interpolator]
    type = PINSFVRhieChowInterpolator
    block = 'reactor pipe pump mixing-plate'
  []
[]

[Variables]
  [pressure]
    type = INSFVPressureVariable
    initial_from_file_var = pressure
    block = 'reactor pipe pump mixing-plate'
  []
  [superficial_vel_x]
    type = PINSFVSuperficialVelocityVariable
    initial_from_file_var = superficial_vel_x
    block = 'reactor pipe pump mixing-plate'
  []
  [superficial_vel_y]
    type = PINSFVSuperficialVelocityVariable
    initial_from_file_var = superficial_vel_y
    block = 'reactor pipe pump mixing-plate'
  []
  [superficial_vel_z]
    type = PINSFVSuperficialVelocityVariable
    initial_from_file_var = superficial_vel_z
    block = 'reactor pipe pump mixing-plate'
  []
  [T]
    type = INSFVEnergyVariable
    initial_from_file_var = T
    block = 'reactor pipe pump mixing-plate'
  []
  [T_ref]
    type = INSFVEnergyVariable
    block = 'reflector'
    initial_from_file_var = T_ref
  []
[]

[AuxVariables]
  [fl]
    type = MooseVariableFVReal
    initial_condition = 1.0
    block = 'reactor pipe pump mixing-plate'
    #initial_from_file_var = fl
  []
  [density_mix]
    type = MooseVariableFVReal
    block = 'reactor pipe pump mixing-plate'
    #initial_from_file_var = density
  []
  [th_cond_mix]
    type = MooseVariableFVReal
    block = 'reactor pipe pump mixing-plate'
    #initial_from_file_var = th_cond
  []
  [cp_mix]
    type = MooseVariableFVReal
    block = 'reactor pipe pump mixing-plate'
    #initial_from_file_var = cp_var
  []
  [darcy_coef]
    type = MooseVariableFVReal
    block = 'reactor pipe pump mixing-plate'
    #initial_from_file_var = darcy_coef
  []
  [fch_coef]
    type = MooseVariableFVReal
    block = 'reactor pipe pump mixing-plate'
    #initial_from_file_var = fch_coef
  []
  [h_DeltaT]
    type = MooseVariableFVReal
    block = 'reactor pipe pump mixing-plate'
  []
  [h_DeltaT_rad_aux]
    type = MooseVariableFVReal
    block = 'reflector'
  []
  [h_DeltaT_rad]
    type = MooseVariableFVReal
    block = 'reflector'
  []
[]

[AuxKernels]
  [compute_fl]
    type = NSLiquidFractionAux
    variable = fl
    temperature = T
    T_liquidus = '${T_liquidus}'
    T_solidus = '${T_solidus}'
    execute_on = 'TIMESTEP_END'
  []
  [rho_out]
    type = ADFunctorElementalAux
    functor = 'rho_mixture'
    variable = 'density_mix'
  []
  [th_cond_out]
    type = ADFunctorElementalAux
    functor = 'k_mixture'
    variable = 'th_cond_mix'
  []
  [cp_out]
    type = ADFunctorElementalAux
    functor = 'cp_mixture'
    variable = 'cp_mix'
  []
  [darcy_out]
    type = ADFunctorElementalAux
    functor = 'Darcy_coefficient'
    variable = 'darcy_coef'
  []
  [fch_out]
    type = ADFunctorElementalAux
    functor = 'Forchheimer_coefficient'
    variable = 'fch_coef'
  []
  [h_DeltaT_out]
    type = ParsedAux
    variable = h_DeltaT
    coupled_variables = 'T'
    expression = '3.*(T-300.)'
  []
  [h_DeltaT_rad_out_pre]
    type = ParsedAux
    variable = h_DeltaT_rad_aux
    coupled_variables = 'T_ref'
    expression = 'T_ref-350.'
  []
  [h_DeltaT_rad_out]
    type = FunctorElementalAux
    functor = 'h_DeltaT_rad_aux'
    variable = h_DeltaT_rad
    factor = 'htc_rad_ref'
  []
[]

[FVKernels]

  ####### MASS EQUATION #######

[mass]
  type = PINSFVMassAdvection
  variable = pressure
  advected_interp_method = 'skewness-corrected'
  velocity_interp_method = 'rc'
  block = 'reactor pipe pump mixing-plate'
[]
####### X-MOMENTUM EQUATION #######

[u_time]
  type = PINSFVMomentumTimeDerivative
  variable = superficial_vel_x
  momentum_component = 'x'
  block = 'reactor pipe pump mixing-plate'
  extra_matrix_tags = 'mass'
[]
[u_advection]
  type = PINSFVMomentumAdvection
  variable = superficial_vel_x
  advected_interp_method = ${advected_interp_method}
  velocity_interp_method = ${velocity_interp_method}
  momentum_component = 'x'
  block = 'reactor pipe pump mixing-plate'
[]
[u_viscosity]
  type = PINSFVMomentumDiffusion
  variable = superficial_vel_x
  momentum_component = 'x'
  block = 'reactor pipe pump mixing-plate'
[]
[u_viscosity_rans]
  type = INSFVMixingLengthReynoldsStress
  variable = superficial_vel_x
  mixing_length = '${fparse 0.07*0.1}'
  momentum_component = 'x'
  block = 'reactor pipe pump mixing-plate'
[]
[u_pressure]
  type = PINSFVMomentumPressure
  variable = superficial_vel_x
  momentum_component = 'x'
  pressure = pressure
  block = 'reactor pipe pump mixing-plate'
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
  Darcy_name = 'DFC_plate'
  Forchheimer_name = 'FFC_plate'
  block = 'mixing-plate'
  consistent_scaling = 100.
[]
#[u_buoyancy]
#  type = PINSFVMomentumBoussinesq
#  variable = superficial_vel_x
#  T_fluid = T
#  ref_temperature = ${T_hx}
#  momentum_component = 'x'
#  #block = 'reactor'
#  alpha_name = ${alpha_b}
#[]

####### Y-MOMENTUM EQUATION #######

[v_time]
 type = PINSFVMomentumTimeDerivative
  variable = superficial_vel_y
  momentum_component = 'y'
  block = 'reactor pipe pump mixing-plate'
  extra_matrix_tags = 'mass'
[]
[v_advection]
  type = PINSFVMomentumAdvection
  variable = superficial_vel_y
  advected_interp_method = ${advected_interp_method}
  velocity_interp_method = ${velocity_interp_method}
  momentum_component = 'y'
  block = 'reactor pipe pump mixing-plate'
[]
[v_viscosity]
  type = PINSFVMomentumDiffusion
  variable = superficial_vel_y
  momentum_component = 'y'
  block = 'reactor pipe pump mixing-plate'
[]
[v_viscosity_rans]
  type = INSFVMixingLengthReynoldsStress
  variable = superficial_vel_y
  mixing_length = '${fparse 0.07*0.1}'
  momentum_component = 'y'
  block = 'reactor pipe pump mixing-plate'
[]
[v_pressure]
  type = PINSFVMomentumPressure
  variable = superficial_vel_y
  momentum_component = 'y'
  pressure = pressure
  block = 'reactor pipe pump mixing-plate'
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

####### Z-MOMENTUM EQUATION #######

[w_time]
  type = PINSFVMomentumTimeDerivative
  variable = superficial_vel_z
  momentum_component = 'z'
  block = 'reactor pipe pump mixing-plate'
  extra_matrix_tags = 'mass'
[]
[w_advection]
  type = PINSFVMomentumAdvection
  variable = superficial_vel_z
  advected_interp_method = ${advected_interp_method}
  velocity_interp_method = ${velocity_interp_method}
  momentum_component = 'z'
  block = 'reactor pipe pump mixing-plate'
[]
[w_viscosity]
  type = PINSFVMomentumDiffusion
  variable = superficial_vel_z
  momentum_component = 'z'
  block = 'reactor pipe pump mixing-plate'
[]
[w_viscosity_rans]
  type = INSFVMixingLengthReynoldsStress
  variable = superficial_vel_z
  mixing_length = '${fparse 0.07*0.1}'
  momentum_component = 'z'
  block = 'reactor pipe pump mixing-plate'
[]
[w_pressure]
  type = PINSFVMomentumPressure
  variable = superficial_vel_z
  momentum_component = 'z'
  pressure = pressure
  block = 'reactor pipe pump mixing-plate'
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

  ####### FUEL ENERGY EQUATION #######

  [heat_time]
    type = PINSFVEnergyTimeDerivative
    variable = T
    is_solid = false
    extra_matrix_tags = 'mass'
  []
  [heat_advection]
    type = PINSFVEnergyAdvection
    variable = T
  []
  [heat_diffusion]
    type = PINSFVEnergyDiffusion
    variable = T
    k = k_mixture
  []
  # [heat_src]
  #   type = FVCoupledForce
  #   variable = T
  #   v = power_density
  #   coef = '${fparse 1.0}'
  #   block = 'reactor'
  # []
  [heat_src]
    type = FVBodyForce
    variable = T
    function = cosine_guess
    value = '${fparse power/0.21757}'  #Volume integral of cosine shape is 0.21757
    block = 'reactor'
  []
  #[heat_sink]
  #  type = PINSFVEnergyAmbientConvection
  #  variable = T
  #  T_fluid = T
  #  T_solid = ${T_hx}
  #  h_solid_fluid = 1e6
  #  is_solid = false
  #  block = 'pump pipe'
  #[]
  [heat_latent_source]
    type = NSFVPhaseChangeSource
    variable = T
    L = ${L}
    liquid_fraction = fl
    T_liquidus = ${T_liquidus}
    T_solidus = ${T_solidus}
  []

  ####### REFLECTOR ENERGY EQUATION #######

  [heat_time_ref]
    type = INSFVEnergyTimeDerivative
    variable = T_ref
    cp = ${cp_ref}
    rho = ${rho_ref}
    extra_matrix_tags = 'mass'
  []
  [heat_diffusion_ref]
    type = FVDiffusion
    variable = T_ref
    coeff = ${k_ref}
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
  [heat-losses-outshield]
    type = FVFunctorConvectiveHeatFluxBC
    variable = T
    T_bulk = T
    T_solid = 300.
    is_solid = false
    heat_transfer_coefficient = 3.
    boundary = 'heat-loss-section-outshield'
  []
  [heated-outshield-pipe]
    type = FVFunctorNeumannBC
    variable = T
    boundary = 'heat-loss-section-outshield'
    functor = heat-input-pipe
  []
  [heat-losses-reflector]
    type = FVFunctorConvectiveHeatFluxBC
    variable = T_ref
    T_bulk = 350.
    T_solid = T_ref
    is_solid = true
    heat_transfer_coefficient = htc_rad_ref
    boundary = 'wall-reflector'
  []
  [heated-reflector-walls]
    type = FVFunctorNeumannBC
    variable = T_ref
    boundary = 'heated-reflector-walls'
    functor = heat-input-ref
  []
[]

[FVInterfaceKernels]
  [convection]
    type = FVConvectionCorrelationInterface
    variable1 = T
    variable2 = T_ref
    boundary = 'wall-reactor-reflector'
    h = htc
    T_solid = T_ref
    T_fluid = T
    subdomain1 = reactor
    subdomain2 = reflector
    wall_cell_is_bulk = true
  []
[]

################################################################################
# MATERIALS
################################################################################

[Functions]
  [heat-input-ref]
    type = ParsedFunction
    expression = '3500.'
  []
  [heat-input-pipe]
    type = ParsedFunction
    expression = '2000.'
  []
  [Re_reactor]
    type = ParsedFunction
    expression = 'flow_hx_bot/Reactor_Area * (2*0.2575) / mu '
    symbol_names = 'flow_hx_bot mu Reactor_Area'
    #symbol_values = 'flow_hx_bot ${mu} ${Reactor_Area}'
    symbol_values = '25.2 ${mu} ${Reactor_Area}'
  []
  [htc]
    type = ParsedFunction
    expression = 'k_liquid/0.2575* 0.023 * Re_reactor^0.8 * Pr^0.3'
    symbol_names = 'k_liquid Re_reactor Pr'
    symbol_values = '${k_liquid} Re_reactor ${Pr}'
  []
  [htc_rad_ref]
    type = ParsedFunction
    expression = '(T_ref_external_wall^2+350.^2)*(T_ref_external_wall+350.)*5.67e-8 / (1/0.18+1/0.35-1.)'   #https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node136.html #e_mgo=0.18 #e_steel=0.35
    symbol_names = 'T_ref_external_wall'
    symbol_values ='T_ref_external_wall'
  []
  [power-density-func]
    type = ParsedFunction
    expression = '${power}/0.21757 * max(0, cos(x*pi/2/1.5))*max(0, cos(y*pi/2/1.5))*max(0, cos(z*pi/2/1.5))'
  []
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
  [cosine_guess]
    type = ParsedFunction
    value = 'max(0, cos(x*pi/2/1.5))*max(0, cos(y*pi/2/1.5))*max(0, cos(z*pi/2/1.5))'
  []
[]

[Materials]
  [generic]
    type = ADGenericFunctorMaterial #defines mu artificially for numerical convergence
    prop_names = 'mu porosity' #it converges to the real mu eventually.
    prop_values = '${mu} ${porosity}'
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
    type = ADGenericVectorFunctorMaterial #defines mu artificially for numerical convergence
    prop_names = 'DFC FFC' #it converges to the real mu eventually.
    prop_values = '${fparse 1*friction} ${fparse 1*friction} ${fparse 1*friction}
                   ${fparse 1*friction} ${fparse 1*friction} ${fparse 1*friction}'
  []
  [friction_material_mixing_plate]
    type = ADGenericVectorFunctorMaterial #defines mu artificially for numerical convergence
    prop_names = 'DFC_plate FFC_plate' #it converges to the real mu eventually.
    prop_values = '${fparse 1*friction} ${fparse 1*friction} ${fparse 1*friction}
                   ${fparse 1*friction} ${fparse 1*friction} ${fparse 1*friction}'
  []
  [ins_fv]
    type = INSFVEnthalpyMaterial
    rho = rho_mixture
    cp = cp_mixture
    temperature = 'T'
  []
  [eff_cp]
    type = NSFVMixtureMaterial
    phase_2_names = '${cp_solid} ${k_solid} ${rho_solid}'
    phase_1_names = '${cp_liquid} ${k_liquid} ${rho_liquid}'
    prop_names = 'cp_mixture k_mixture rho_mixture'
    phase_1_fraction = fl
  []
  [mushy_zone_resistance]
    type = INSFVMushyPorousFrictionMaterial
    liquid_fraction = 'fl'
    mu = '${mu}'
    rho_l = '${rho_liquid}'
    dendrite_spacing_scaling = 1e-2
  []
[]

[Problem]
  error_on_jacobian_nonzero_reallocation = true
  extra_tag_matrices = 'mass'
  type = NavierStokesProblem
  mass_matrix = 'mass'
  schur_fs_index = '1'
[]

################################################################################
# EXECUTION / SOLVE
################################################################################
[Preconditioning]
  active = FSP
  [FSP]
    type = FSP
    topsplit = 'top'
    [top]
      splitting = 'up temperatures'
      splitting_type = multiplicative
      petsc_options_iname = '-ksp_type'
      petsc_options_value = 'preonly'
    []
      [up]
        splitting = 'us p'
        splitting_type  = schur
        petsc_options_iname = '-pc_fieldsplit_schur_fact_type  -pc_fieldsplit_schur_precondition -ksp_gmres_restart -ksp_rtol -ksp_type -ksp_atol'
        petsc_options_value = 'full                            selfp                             300                1e-5      fgmres 1e-9'
        vars = 'superficial_vel_x superficial_vel_y superficial_vel_z pressure'
      []
        [us]
          vars = 'superficial_vel_x superficial_vel_y superficial_vel_z'
          splitting = 'u v w'
          splitting_type = symmetric_multiplicative
          petsc_options = '-ksp_monitor'
          petsc_options_iname = ' -ksp_gmres_restart -ksp_rtol -ksp_type -ksp_atol'
          petsc_options_value = '300                1e-5      fgmres 1e-9'
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
          petsc_options = '-ksp_converged_reason -ksp_monitor'
          petsc_options_iname = '-ksp_type -ksp_gmres_restart -ksp_rtol -pc_type -ksp_pc_side'
          petsc_options_value = 'fgmres     300                1e-5     lu       right'
        []
      [temperatures]
        splitting_type = additive
        splitting = 'T T_ref'
        petsc_options_iname = ' -ksp_gmres_restart -ksp_rtol -ksp_type -ksp_atol'
        petsc_options_value = '300                1e-5      fgmres 1e-9'
        vars = 'T T_ref'
      []
        [T]
          vars = 'T'
          petsc_options = '-ksp_monitor_true_residual -ksp_converged_reason'
          petsc_options_iname = '-pc_type -ksp_pc_side -ksp_type -ksp_rtol'
          petsc_options_value = 'lu       right        gmres     1e-5'
        []
        [T_ref]
          vars = 'T_ref'
          petsc_options = '-ksp_monitor_true_residual -ksp_converged_reason'
          petsc_options_iname = '-pc_type -ksp_pc_side -ksp_type -ksp_rtol'
          petsc_options_value = 'lu       right        gmres     1e-5'
        []
  []
[]

[Executioner]
  type = Transient

  # Time-stepping parameters
  start_time = 0.0
  end_time = 864000

  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 10
    dt = 1e-3
    timestep_limiting_postprocessor = 'dt_limit'
  []

  # Solver parameters
  solve_type = 'NEWTON'
  line_search = 'none'
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-3
  nl_max_its = 10 # fail early and try again with a shorter time step
  l_max_its = 80
  # automatic_scaling = true
[]

[Debug]
  show_var_residual_norms = true
[]

################################################################################
# SIMULATION OUTPUTS
################################################################################

[Outputs]
  csv = true
  [restart]
    type = Exodus
    execute_on = 'timestep_end final'
  []
  # Reduce base output
  print_linear_converged_reason = true
  print_linear_residuals = false
  print_nonlinear_converged_reason = true
  hide = 'T_hx_top pdrop dt_limit'
[]

[Postprocessors]
  #[max_v]
  #  type = ElementExtremeValue
  #  variable = superficial_vel_x
  #  value_type = max
  #  block = 'reactor pump pipe'
  #[]
  #[flow_hx_bot]
  #  type = VolumetricFlowRate
  #  boundary = 'reactor_bot'
  #  vel_x = superficial_vel_x
  #  vel_y = superficial_vel_y
  #  advected_quantity = ${rho}
  #[]
  #[flow_hx_top]
  #  type = VolumetricFlowRate
  #  boundary = 'reactor_top'
  #  vel_x = superficial_vel_x
  #  vel_y = superficial_vel_y
  #  advected_quantity = ${rho}
  #[]
  [T_hx_bot]
    type = SideAverageValue
    boundary = 'reactor_bot'
    variable = T
  []
  [T_hx_top]
    type = SideAverageValue
    boundary = 'reactor_top'
    variable = T
  []
  [T_cooled_inlet]
    type = SideAverageValue
    boundary = 'inlet-cooled-wall-pipes'
    variable = T
  []
  [T_cooled_outlet]
    type = SideAverageValue
    boundary = 'outlet-cooled-wall-pipes'
    variable = T
  []
  [T_ref_external_wall]
    type = SideAverageValue
    boundary = 'wall-reflector'
    variable = T_ref
  []
  [pdrop]
    type = PressureDrop
    pressure = pressure
    upstream_boundary = 'reactor_bot'
    downstream_boundary = 'mixing-plate-downstream'
    boundary = 'mixing-plate-downstream reactor_bot'
  []
  [heat-loss-pipe]
    type = SideIntegralVariablePostprocessor
    boundary = heat-loss-section-outshield
    variable = h_DeltaT
  []
  [heat-loss-reflector]
    type = SideIntegralVariablePostprocessor
    boundary = wall-reflector
    variable = h_DeltaT_rad
  []
  [heat-input-ref]
    type = FunctionSideIntegral
    boundary = heated-reflector-walls
    function = heat-input-ref
  []
  [heat-input-pipe]
    type = FunctionSideIntegral
    boundary = heat-loss-section-outshield
    function = heat-input-pipe
  []
  [reactor-power]
    type = ElementIntegralFunctorPostprocessor
    block = 'reactor'
    functor = 'power-density-func'
  []
  #[area]
  #  type = AreaPostprocessor
  #  boundary = 'wall-reflector'
  #[]
  [dt_limit]
    type = Receiver
    default = 10000
  []
[]
