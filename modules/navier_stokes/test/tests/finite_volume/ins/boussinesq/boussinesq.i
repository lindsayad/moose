mu = 1
rho = 1
k = 1
cp = 1
alpha = 1
advected_interp_method = 'upwind'
# rayleigh=1e3
# hot_temp=${rayleigh}
hot_temp=1
temp_ref=${fparse hot_temp / 2.}
l=1e2

[GlobalParams]
  rhie_chow_user_object = 'rc'
  velocity_interp_method = 'rc'
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = u
    v = v
    pressure = pressure
    # mu = 'mu'
    # rho = 'rho'
    # gravity = '0 -1 0'
    # alpha_b = ${alpha}
    # T_fluid = 'T'
  []
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${l}
    ymin = 0
    ymax = ${l}
    nx = 8
    ny = 8
  []
[]

[Variables]
  [u]
    type = INSFVVelocityVariable
    initial_condition = 1e-12
  []
  [v]
    type = INSFVVelocityVariable
    initial_condition = 1e-12
  []
  [pressure]
    type = INSFVPressureVariable
  []
  [T]
    type = INSFVEnergyVariable
    # scaling = 1e-4
  []
  [lambda]
    family = SCALAR
    order = FIRST
  []
[]

[AuxVariables]
  [U]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
  [vel_x]
    order = FIRST
    family = MONOMIAL
  []
  [vel_y]
    order = FIRST
    family = MONOMIAL
  []
  [viz_T]
    order = FIRST
    family = MONOMIAL
  []
[]

[AuxKernels]
  [mag]
    type = VectorMagnitudeAux
    variable = U
    x = u
    y = v
    execute_on = 'initial timestep_end'
  []
  [vel_x]
    type = ParsedAux
    variable = vel_x
    function = 'u'
    execute_on = 'initial timestep_end'
    args = 'u'
  []
  [vel_y]
    type = ParsedAux
    variable = vel_y
    function = 'v'
    execute_on = 'initial timestep_end'
    args = 'v'
  []
  [viz_T]
    type = ParsedAux
    variable = viz_T
    function = 'T'
    execute_on = 'initial timestep_end'
    args = 'T'
  []
[]

[FVKernels]
  [mass]
    type = INSFVMassAdvection
    variable = pressure
    advected_interp_method = ${advected_interp_method}
    rho = ${rho}
  []
  [mean_zero_pressure]
    type = FVScalarLagrangeMultiplier
    variable = pressure
    lambda = lambda
  []

  [u_time]
    type = INSFVMomentumTimeDerivative
    rho = ${rho}
    momentum_component = 'x'
    variable = u
  []
  [u_advection]
    type = INSFVMomentumAdvection
    variable = u
    advected_interp_method = ${advected_interp_method}
    rho = ${rho}
    momentum_component = 'x'
    # velocity_interp_method = 'average'
  []
  [u_viscosity]
    type = INSFVMomentumDiffusion
    variable = u
    mu = ${mu}
    momentum_component = 'x'
  []
  [u_pressure]
    type = INSFVMomentumPressure
    variable = u
    momentum_component = 'x'
    pressure = pressure
  []
  [u_buoyancy]
    type = INSFVMomentumBoussinesq
    variable = u
    T_fluid = T
    gravity = '0 -1 0'
    rho = ${rho}
    ref_temperature = ${temp_ref}
    momentum_component = 'x'
  []
  [u_gravity]
    type = INSFVMomentumGravity
    variable = u
    gravity = '0 -1 0'
    rho = ${rho}
    momentum_component = 'x'
  []

  [v_time]
    type = INSFVMomentumTimeDerivative
    rho = ${rho}
    momentum_component = 'y'
    variable = v
  []
  [v_advection]
    type = INSFVMomentumAdvection
    variable = v
    advected_interp_method = ${advected_interp_method}
    rho = ${rho}
    momentum_component = 'y'
    # velocity_interp_method = 'average'
  []
  [v_viscosity]
    type = INSFVMomentumDiffusion
    variable = v
    mu = ${mu}
    momentum_component = 'y'
  []
  [v_pressure]
    type = INSFVMomentumPressure
    variable = v
    momentum_component = 'y'
    pressure = pressure
  []
  [v_buoyancy]
    type = INSFVMomentumBoussinesq
    variable = v
    T_fluid = T
    gravity = '0 -1 0'
    rho = ${rho}
    ref_temperature = ${temp_ref}
    momentum_component = 'y'
  []
  [v_gravity]
    type = INSFVMomentumGravity
    variable = v
    gravity = '0 -1 0'
    rho = ${rho}
    momentum_component = 'y'
  []

  [temp_time]
    type = INSFVEnergyTimeDerivative
    variable = T
    rho = ${rho}
  []
  [temp_conduction]
    type = FVDiffusion
    coeff = 'k'
    variable = T
  []
  [temp_advection]
    type = INSFVEnergyAdvection
    variable = T
    advected_interp_method = ${advected_interp_method}
    # velocity_interp_method = 'average'
  []
[]

[FVBCs]
  [top_x]
    type = INSFVNoSlipWallBC
    variable = u
    boundary = 'top'
    function = 'lid_function'
  []

  [no_slip_x]
    type = INSFVNoSlipWallBC
    variable = u
    boundary = 'left right bottom'
    function = 0
  []

  [no_slip_y]
    type = INSFVNoSlipWallBC
    variable = v
    boundary = 'left right top bottom'
    function = 0
  []

  [T_hot]
    type = FVDirichletBC
    variable = T
    boundary = left
    value = ${hot_temp}
  []

  [T_cold]
    type = FVDirichletBC
    variable = T
    boundary = right
    value = 0
  []
[]

[Problem]
  error_on_jacobian_nonzero_reallocation = false
[]

[Materials]
  [const_functor]
    type = ADGenericFunctorMaterial
    prop_names = 'alpha_b cp k mu rho'
    prop_values = '${alpha} ${cp} ${k} ${mu} ${rho}'
  []
  [ins_fv]
    type = INSFVEnthalpyMaterial
    temperature = 'T'
    rho = ${rho}
  []
[]

[Functions]
  [lid_function]
    type = ParsedFunction
    value = '4*x*(1-x)'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  # type = Steady
  type = Transient
  solve_type = 'NEWTON'
  petsc_options = '-pc_svd_monitor'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'svd'
  # petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
  # petsc_options_value = 'asm      300                lu           NONZERO'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-7
  nl_max_its = 10
  line_search = 'none'
  l_max_its = 20
  end_time = 1e9
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    dt = 1e-3
    growth_factor = 1.2
  []
[]

[Outputs]
  exodus = true
[]

# [Postprocessors]
#   [Rayleigh]
#     type = RayleighNumber
#     beta = ${alpha}
#     T_hot = ${hot_temp}
#     T_cold = 0
#     rho_ave = ${rho}
#     l = ${l}
#     mu_ave = ${mu}
#     k_ave = ${k}
#     cp_ave = ${cp}
#     gravity_magnitude = 1
#   []
# []
