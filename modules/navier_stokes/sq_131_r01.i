# Pronghorn assessment
# Square Enclosure Test 131
#
hot_temp = 336.84
cold_temp = 297.13
#
mu = 1
rho = 1
k = 1
cp = 1
alpha = 1
vel = 'velocity'
velocity_interp_method = 'rc'
advected_interp_method = 'upwind'
#rayleigh=1e3
#hot_temp=${rayleigh}
#temp_ref=${fparse hot_temp / 2.}
temp_ref = 316.99
#
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 0.0635
    ymin = 0
    ymax = 0.0635
    nx = 100
    ny = 100
  []
[]

[Variables]
  [u]
    type = INSFVVelocityVariable
  []
  [v]
    type = INSFVVelocityVariable
  []
  [pressure]
    type = INSFVPressureVariable
  []
  [T]
    type = INSFVEnergyVariable
    scaling = 1e-4
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
  [./flux_x]
      order = FIRST
      family = MONOMIAL
  [../]
  [./flux_y]
      order = FIRST
      family = MONOMIAL
  [../]
  [./QLeft]
      order = FIRST
      family = MONOMIAL
  [../]
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
    vel = ${vel}
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    u = u
    v = v
    pressure = pressure
    mu = ${mu}
    rho = ${rho}
  []
  [mean_zero_pressure]
    type = FVScalarLagrangeMultiplier
    variable = pressure
    lambda = lambda
  []

  [u_time]
    type = INSFVMomentumTimeDerivative
    variable = 'u'
    rho = ${rho}
  []
  [u_advection]
    type = INSFVMomentumAdvection
    variable = u
    advected_quantity = 'rhou'
    vel = ${vel}
    velocity_interp_method = ${velocity_interp_method}
    advected_interp_method = ${advected_interp_method}
    pressure = pressure
    u = u
    v = v
    mu = ${mu}
    rho = ${rho}
  []
  [u_viscosity]
    type = FVDiffusion
    variable = u
    coeff = ${mu}
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
    variable = v
    rho = ${rho}
  []
  [v_advection]
    type = INSFVMomentumAdvection
    variable = v
    advected_quantity = 'rhov'
    vel = ${vel}
    velocity_interp_method = ${velocity_interp_method}
    advected_interp_method = ${advected_interp_method}
    pressure = pressure
    u = u
    v = v
    mu = ${mu}
    rho = ${rho}
  []
  [v_viscosity]
    type = FVDiffusion
    variable = v
    coeff = ${mu}
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
    cp_name = 'cp'
  []
  [temp_conduction]
    type = FVDiffusion
    coeff = 'k'
    variable = T
  []
  [temp_advection]
    type = INSFVEnergyAdvection
    variable = T
    vel = ${vel}
    velocity_interp_method = ${velocity_interp_method}
    advected_interp_method = ${advected_interp_method}
    pressure = pressure
    u = u
    v = v
    mu = ${mu}
    rho = ${rho}
  []
[]

[FVBCs]
##  [top_x]
##    type = INSFVNoSlipWallBC
##    variable = u
##    boundary = 'left right top bottom'
##    value = 0
##  []
#
  [no_slip_x]
    type = INSFVNoSlipWallBC
##    type = FVDirichletBC
    variable = u
    boundary = 'left right top bottom'
    function = 0
  []
#
  [no_slip_y]
    type = INSFVNoSlipWallBC
##    type = FVDirichletBC
    variable = v
    boundary = 'left right top bottom'
    function = 0
  []
#
  [T_hot]
    type = FVDirichletBC
    variable = T
    boundary = left
    value = ${hot_temp}
  []
#
  [T_cold]
    type = FVDirichletBC
    variable = T
    boundary = right
    value = ${cold_temp}
  []
[]
#
[Materials]
  [const_functor]
    type = ADGenericFunctorMaterial
    prop_names = 'alpha_b cp k'
    prop_values = '${alpha} ${cp} ${k}'
  []
  [ins_fv]
    type = INSFVMaterial
    u = 'u'
    v = 'v'
    pressure = 'pressure'
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
#[Preconditioning]
#  [./SMP_PJFNK]
#    type = SMP
#    full = true
#    solve_type = 'PJFNK'
#  [../]
#[]
#
[Executioner]
  type = Transient
  dtmin = 0.00001
  start_time = 0.0
  end_time = 100.0
  num_steps = 25000
  steady_state_detection = true
  steady_state_tolerance = 1e-6
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    growth_factor = 1.5
    dt = 0.01
  []
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
[]

[Outputs]
print_linear_residuals = false
perf_graph = true
interval = 1
execute_on = 'initial timestep_end'
#exodus = true
#  output_material_properties = true
[./out]
  type = Exodus
  output_material_properties = true
[../]
##[./checkpoint]
##  type = Checkpoint
##  num_files = 1
##[../]
[./csv]
  type = CSV
  interval = 1
[../]
[./console]
  type = Console
  output_linear = true
[../]
[]
#
[Postprocessors]
  [./timestep]
    type = TimestepSize
  [../]
  [./timemax]
    type = TimeExtremeValue
    postprocessor = timestep
    value_type = max
  [../]
  [./timemin]
    type = TimeExtremeValue
    postprocessor = timestep
    value_type = min
  [../]
  [./max]
    type = ElementExtremeValue
    variable = T
    value_type = max
  [../]
  [./min]
    type = ElementExtremeValue
    variable = T
    value_type = min
  [../]
  [./hotwallflux]
    type = SideFluxAverage
    variable = T
    diffusivity = 0.0265
    boundary = 'left'
  [../]
  [./hotwallintflux]
    type = SideFluxIntegral
    variable = T
    diffusivity = 0.0265
    boundary = 'left'
  [../]
  [./coldwallflux]
    type = SideFluxAverage
    variable = T
    diffusivity = 0.0265
    boundary = 'right'
  [../]
  [./coldwallintflux]
    type = SideFluxIntegral
    variable = T
    diffusivity = 0.0265
    boundary = 'right'
  [../]
  [./topwallflux]
    type = SideFluxAverage
    variable = T
    diffusivity = 0.0265
    boundary = 'top'
  [../]
  [./bottomwallflux]
    type = SideFluxAverage
    variable = T
    diffusivity = 0.0265
    boundary = 'bottom'
  [../]
  [./Tave]
    type = ElementAverageValue
    variable = T
  [../]
[]
#
#
[VectorPostprocessors]
  [./fx_left]
    type = SideValueSampler
    variable = flux_x
    boundary = 'left'
    sort_by  = y
  [../]
  [./fx_right]
    type = SideValueSampler
    variable = flux_x
    boundary = 'right'
    sort_by  = y
  [../]
  [./T_mid]
    type = LineValueSampler
    variable = T
    sort_by  = y
    start_point = '0.03175 0 0'
    end_point   = '0.03175 0.0635 0'
    num_points  = 21
  [../]
  [./QLeftVPP]
    type = SideValueSampler
    variable = QLeft
    boundary = 'left'
    sort_by  = y
  [../]
[]
