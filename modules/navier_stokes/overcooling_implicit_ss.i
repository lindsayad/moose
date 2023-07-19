################################################################################
## Molten Salt Fast Reactor - Euratom EVOL + Rosatom MARS Design              ##
## Pronghorn input file to initialize velocity fields                         ##
## This runs a slow relaxation to steady state while ramping down the fluid   ##
## viscosity.                                                                 ##
################################################################################

# Material properties fuel
rho = 3279. # density [kg / m^3]  (@1000K)
mu = 0.005926 # viscosity [Pa s]
rho_solid = 3279.7
rho_liquid = 3279.7
k_solid = 0.38
k_liquid = 0.38
cp_solid = 640.
cp_liquid = 640.
L = 4e5
T_liquidus = 798. #812.
T_solidus =  790. #785.

# Material properties reflector
k_ref = 30.
cp_ref = 880.
rho_ref = 3580.

power = 25e3

# Mass flow rate tuning
friction = 11.0 # [kg / m^4]
pump_force = '${fparse 0.25*4.0e6}' # [N / m^3]
porosity = 1.0
Reactor_Area = ${fparse 3.14159*0.2575*0.2575}
Pr = ${fparse mu*cp_liquid/k_solid}

# Numerical scheme parameters
advected_interp_method = 'upwind'
velocity_interp_method = 'average'

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

  #temperature = T
[]

################################################################################
# GEOMETRY
################################################################################
[Mesh]
  [fmg]
    type = FileMeshGenerator
    use_for_exodus_restart = true
    file = overcooling_implicit_ss_rc_restart.e #overcooling_implicit_ss_2.e overcooling_implicit_ss_rc_restart.e
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

[GlobalParams]
    initial_from_file_timestep = LATEST
[]

[Variables]
  [pressure]
    type = INSFVPressureVariable
    initial_from_file_var = pressure
    block = 'reactor pipe pump mixing-plate'
    initial_from_file_timestep = LATEST
  []
  [superficial_vel_x]
    type = PINSFVSuperficialVelocityVariable
    initial_from_file_var = superficial_vel_x
    block = 'reactor pipe pump mixing-plate'
    initial_from_file_timestep = LATEST
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
  #[lambda]
  #  family = SCALAR
  #  order = FIRST
  #  block = 'reactor pipe pump mixing-plate'
  #  initial_from_file_var = lambda
  #[]
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
  # [T]
  #   type = MooseVariableFVReal
  #   initial_from_file_var = T
  #   block = 'reactor pipe pump mixing-plate'
  # []
  # [T_ref]
  #   type = MooseVariableFVReal
  #   block = 'reflector'
  #   initial_from_file_var = T_ref
  # []
  [a_u_var]
    type = MooseVariableFVReal
    initial_from_file_var = a_u_var
  []
  [a_v_var]
    type = MooseVariableFVReal
    initial_from_file_var = a_v_var
  []
  [a_w_var]
    type = MooseVariableFVReal
    initial_from_file_var = a_w_var
  []
  [T_debug]
    type = MooseVariableFVReal
    block = 'reactor pipe pump mixing-plate'
  []
  [fl]
    type = MooseVariableFVReal
    #initial_condition = 1.0
    block = 'reactor pipe pump mixing-plate'
    initial_from_file_var = fl
  []
  [density_mix]
    type = MooseVariableFVReal
    block = 'reactor pipe pump mixing-plate'
    initial_from_file_var = density_mix
  []
  [th_cond_mix]
    type = MooseVariableFVReal
    block = 'reactor pipe pump mixing-plate'
    initial_from_file_var = th_cond_mix
  []
  [cp_mix]
    type = MooseVariableFVReal
    block = 'reactor pipe pump mixing-plate'
    initial_from_file_var = cp_mix
  []
  [darcy_coef]
    type = MooseVariableFVReal
    block = 'reactor pipe pump mixing-plate'
    initial_from_file_var = darcy_coef
  []
  [fch_coef]
    type = MooseVariableFVReal
    block = 'reactor pipe pump mixing-plate'
    initial_from_file_var = fch_coef
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
  # [comp_a_u]
  #   type = FunctorElementalAux
  #   functor = 'ax'
  #   variable = 'a_u_var'
  #   block = 'reactor pipe pump mixing-plate'
  #   execute_on = 'timestep_end'
  # []
  # [comp_a_v]
  #   type = FunctorElementalAux
  #   functor = 'ay'
  #   variable = 'a_v_var'
  #   block = 'reactor pipe pump mixing-plate'
  #   execute_on = 'timestep_end'
  # []
  # [comp_a_w]
  #   type = FunctorElementalAux
  #   functor = 'az'
  #   variable = 'a_w_var'
  #   block = 'reactor pipe pump mixing-plate'
  #   execute_on = 'timestep_end'
  # []
  [compute_T_debug]
    type = DebugResidualAux
    variable = T_debug
    debug_variable = T
    execute_on = NONLINEAR
  []
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
  velocity_interp_method = ${velocity_interp_method}
  block = 'reactor pipe pump mixing-plate'
[]
#[mean_zero_pressure]
#  #type = FVIntegralValueConstraint
#  type = FVPointValueConstraint
#  variable = pressure
#  lambda = lambda
#  point = '0. 0. 0.'
#  block = 'reactor pipe pump mixing-plate'
#[]
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
[u_friction_mushyzone]
  type = INSFVMomentumFriction
  variable = superficial_vel_x
  momentum_component = 'x'
  linear_coef_name = 'Darcy_coefficient'
  quadratic_coef_name = 'Forchheimer_coefficient'
  block = 'reactor pipe'
[]

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
[v_friction_mushyzone]
  type = INSFVMomentumFriction
  variable = superficial_vel_y
  momentum_component = 'y'
  linear_coef_name = 'Darcy_coefficient'
  quadratic_coef_name = 'Forchheimer_coefficient'
  block = 'reactor pipe'
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
[w_friction_mushyzone]
  type = INSFVMomentumFriction
  variable = superficial_vel_z
  momentum_component = 'z'
  linear_coef_name = 'Darcy_coefficient'
  quadratic_coef_name = 'Forchheimer_coefficient'
  block = 'reactor pipe'
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
    block = 'reactor pipe'
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
    heat_transfer_coefficient = 50
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
  [heated-inshield-pipe]
    type = FVNeumannBC
    variable = T
    boundary = 'heat-loss-section-inshield'
    value = 1000.
  []
[]

[FVInterfaceKernels]
  [convection]
    type = FVConvectionCorrelationInterface
    variable1 = T
    variable2 = T_ref
    boundary = 'wall-reactor-full'
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
    expression = '0.'
  []
  [heat-input-pipe]
    type = ParsedFunction
    expression = '0.'
  []
  [Re_reactor]
    type = ParsedFunction
    expression = 'if(flow_hx_bot<1., 1.058525e+04, flow_hx_bot/Reactor_Area * (2*0.2575) / mu)'
    symbol_names = 'flow_hx_bot mu Reactor_Area'
    #symbol_values = 'flow_hx_bot ${mu} ${Reactor_Area}'
    symbol_values = 'flow_hx_bot ${mu} ${Reactor_Area}'
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
    type = ParsedFunction
    expression = mu*(0.1*exp(-3*t)+1)
    symbol_names = 'mu'
    symbol_values = ${mu}
  []
  [mu_x]
    type = ParsedFunction
    expression = 'if(x > -0.08 & x < 0.2 , 10.*mu*(0.1*exp(-3*t)+1),mu*(0.1*exp(-3*t)+1))'
    symbol_names = 'mu'
    symbol_values = ${mu}
  []
  [cosine_guess]
    type = ParsedFunction
    expression = 'max(0, cos(x*pi/2/1.5))*max(0, cos(y*pi/2/1.5))*max(0, cos(z*pi/2/1.5))'
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
  extra_tag_matrices = 'mass L'
  type = NavierStokesProblem
  mass_matrix = 'mass'
  L_matrix = 'L'
  use_mass_matrix_for_scaling = false
  commute_lsc = false
  material_coverage_check = false
  kernel_coverage_check = false
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
        petsc_options = '-ksp_monitor_true_residual'
        petsc_options_iname = '-pc_fieldsplit_schur_fact_type  -pc_fieldsplit_schur_precondition -ksp_gmres_restart -ksp_rtol -ksp_type -ksp_atol'
        petsc_options_value = 'full                            self                             300                 1e-2      fgmres    7e-5'
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
      [temperatures]
        splitting_type = additive
        splitting = 'T T_ref'
        petsc_options_iname = ' -ksp_gmres_restart -ksp_rtol -ksp_type -ksp_atol'
        petsc_options_value = '300                1e-5      fgmres 1e-9'
        vars = 'T T_ref'
        petsc_options = '-ksp_monitor'
      []
        [T]
          vars = 'T'
          petsc_options_iname = '-pc_type -ksp_pc_side -ksp_type -ksp_rtol'
          petsc_options_value = 'lu       right        gmres     1e-5'
        []
        [T_ref]
          vars = 'T_ref'
          petsc_options_iname = '-pc_type -ksp_pc_side -ksp_type -ksp_rtol'
          petsc_options_value = 'lu       right        gmres     1e-5'
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
  end_time = 10000.
  dtmax = 100
  steady_state_detection = true

  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 7
    dt = 1.
    # timestep_limiting_postprocessor = 'dt_limit'
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
  [testsss]
    type = Exodus
    execute_on = 'initial timestep_end final'
  []
  # Reduce base output
  print_linear_converged_reason = false
  print_linear_residuals = false
  print_nonlinear_converged_reason = false
  hide = 'T_hx_top pdrop dt_limit'
  checkpoint = true
[]

[Postprocessors]
  #[max_v]
  #  type = ElementExtremeValue
  #  variable = superficial_vel_x
  #  value_type = max
  #  block = 'reactor pump pipe'
  #[]
  [flow_hx_bot]
    type = VolumetricFlowRate
    boundary = 'reactor_bot'
    vel_x = superficial_vel_x
    vel_y = superficial_vel_y
    advected_quantity = ${rho}
  []
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
    execute_on = 'initial timestep_end'
  []
  [T_hx_top]
    type = SideAverageValue
    boundary = 'reactor_top'
    variable = T
    execute_on = 'initial timestep_end'
  []
  [T_cooled_inlet]
    type = SideAverageValue
    boundary = 'inlet-cooled-wall-pipes'
    variable = T
    execute_on = 'initial timestep_end'
  []
  [T_cooled_outlet]
    type = SideAverageValue
    boundary = 'outlet-cooled-wall-pipes'
    variable = T
    execute_on = 'initial timestep_end'
  []
  [T_ref_external_wall]
    type = SideAverageValue
    boundary = 'wall-reflector'
    variable = T_ref
    execute_on = 'initial timestep_end'
  []
  [pdrop]
    type = PressureDrop
    pressure = pressure
    upstream_boundary = 'reactor_bot'
    downstream_boundary = 'mixing-plate-downstream'
    boundary = 'mixing-plate-downstream reactor_bot'
    execute_on = 'initial timestep_end'
  []
  [heat-loss-pipe]
    type = SideIntegralVariablePostprocessor
    boundary = heat-loss-section-outshield
    variable = h_DeltaT
    execute_on = 'initial timestep_end'
  []
  [heat-loss-reflector]
    type = SideIntegralVariablePostprocessor
    boundary = wall-reflector
    variable = h_DeltaT_rad
    execute_on = 'initial timestep_end'
  []
  [heat-input-ref]
    type = FunctionSideIntegral
    boundary = heated-reflector-walls
    function = heat-input-ref
    execute_on = 'initial timestep_end'
  []
  [heat-input-pipe]
    type = FunctionSideIntegral
    boundary = heat-loss-section-outshield
    function = heat-input-pipe
    execute_on = 'initial timestep_end'
  []
  [reactor-power]
    type = ElementIntegralFunctorPostprocessor
    block = 'reactor'
    functor = 'power-density-func'
    execute_on = 'initial timestep_end'
  []
  #[area]
  #  type = AreaPostprocessor
  #  boundary = 'wall-reflector'
  #[]
  [Re_reactor]
    type = FunctionValuePostprocessor
    function = Re_reactor
    execute_on = 'initial timestep_end'
  []
  [min_T]
    type = ElementExtremeValue
    variable = T
    value_type = min
    block = 'reactor pipe mixing-plate'
    execute_on = timestep_end
  []
  [htc]
    type = FunctionValuePostprocessor
    function = htc
    execute_on = 'initial timestep_end'
  []
  [dt_limit]
    type = Receiver
    default = 100.
  []
[]
