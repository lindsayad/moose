mu = 1.1
rho = 1.1
Vin = 1.0

width = 1		# channel half-width
length = 10

k = FIRST
k_minus_one = CONSTANT
pressure_family = MONOMIAL

alpha = 6

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${length}
    ymin = ${fparse -width}
    ymax = ${width}
    nx = 100
    ny = 20

#     elem_type = TRI6
    elem_type = QUAD9
  []
[]

[Variables]
  [vel_x]
    family = L2_HIERARCHIC
    order = ${k}
  []
  [vel_y]
    family = L2_HIERARCHIC
    order = ${k}
  []
  [pressure]
    family = ${pressure_family}
    order = ${k_minus_one}
  []
  [vel_bar_x]
    family = SIDE_HIERARCHIC
    order = ${k}
  []
  [vel_bar_y]
    family = SIDE_HIERARCHIC
    order = ${k}
  []
  [pressure_bar]
    family = SIDE_HIERARCHIC
    order = ${k}
  []
[]

[HDGKernels]
  [momentum_x_advection]
    type = AdvectionIPHDGKernel
    variable = vel_x
    face_variable = vel_bar_x
    velocity = 'velocity'
    coeff = ${rho}
  []
  [momentum_x_diffusion]
    type = NavierStokesStressIPHDGKernel
    variable = vel_x
    face_variable = vel_bar_x
    diffusivity = 'mu'
    pressure_variable = pressure
    pressure_face_variable = pressure_bar
    component = 0
    alpha = ${alpha}
  []
  [momentum_y_advection]
    type = AdvectionIPHDGKernel
    variable = vel_y
    face_variable = vel_bar_y
    velocity = 'velocity'
    coeff = ${rho}
  []
  [momentum_y_diffusion]
    type = NavierStokesStressIPHDGKernel
    variable = vel_y
    face_variable = vel_bar_y
    diffusivity = 'mu'
    pressure_variable = pressure
    pressure_face_variable = pressure_bar
    component = 1
    alpha = ${alpha}
  []
  [pressure]
    type = MassContinuityIPHDGKernel
    variable = pressure
    face_variable = pressure_bar
    interior_velocity_vars = 'vel_x vel_y'
    face_velocity_functors = 'vel_bar_x vel_bar_y'
  []
[]

[Functions]
  [u_inlet]
    type = ConstantFunction
    value = ${Vin}
  []
[]

[BCs]
  #
  # dirichlet
  #
  [momentum_x_advection_dirichlet]
    type = AdvectionIPHDGDirichletBC
    boundary = left
    face_variable = vel_bar_x
    functor = u_inlet
    variable = vel_x
    velocity = velocity
    coeff = ${rho}
  []
  [momentum_x_diffusion_dirichlet]
    type = NavierStokesStressIPHDGDirichletBC
    boundary = left
    variable = vel_x
    face_variable = vel_bar_x
    pressure_variable = pressure
    pressure_face_variable = pressure_bar
    alpha = ${alpha}
    functor = 'u_inlet'
    diffusivity = 'mu'
    component = 0
  []
  [momentum_y_diffusion_dirichlet]
    type = NavierStokesStressIPHDGDirichletBC
    boundary = left
    variable = vel_y
    face_variable = vel_bar_y
    pressure_variable = pressure
    pressure_face_variable = pressure_bar
    alpha = ${alpha}
    functor = '0'
    diffusivity = 'mu'
    component = 1
  []
  [mass_dirichlet]
    type = MassContinuityIPHDGBC
    face_variable = pressure_bar
    variable = pressure
    boundary = left
    face_velocity_functors = 'u_inlet 0'
    interior_velocity_vars = 'vel_x vel_y'
  []

  #
  # walls
  #
  [momentum_x_diffusion_walls]
    type = NavierStokesStressIPHDGDirichletBC
    boundary = 'top bottom'
    variable = vel_x
    face_variable = vel_bar_x
    pressure_variable = pressure
    pressure_face_variable = pressure_bar
    alpha = ${alpha}
    functor = '0'
    diffusivity = 'mu'
    component = 0
  []
  [momentum_y_diffusion_walls]
    type = NavierStokesStressIPHDGDirichletBC
    boundary = 'top bottom'
    variable = vel_y
    face_variable = vel_bar_y
    pressure_variable = pressure
    pressure_face_variable = pressure_bar
    alpha = ${alpha}
    functor = 0
    diffusivity = 'mu'
    component = 1
  []
  [mass_walls]
    type = MassContinuityIPHDGBC
    face_variable = pressure_bar
    variable = pressure
    boundary = 'top bottom'
    face_velocity_functors = '0 0'
    interior_velocity_vars = 'vel_x vel_y'
  []

  #
  # Neumann
  #
  [momentum_x_advection_neumann]
    type = AdvectionIPHDGOutflowBC
    boundary = 'right'
    constrain_lm = false
    face_variable = vel_bar_x
    variable = vel_x
    velocity = velocity
    coeff = ${rho}
  []
  [momentum_y_advection_neumann]
    type = AdvectionIPHDGOutflowBC
    boundary = 'right'
    constrain_lm = false
    face_variable = vel_bar_y
    variable = vel_y
    velocity = velocity
    coeff = ${rho}
  []
  [momentum_x_diffusion_neumann]
    type = NavierStokesStressIPHDGPrescribedTractionBC
    boundary = 'right'
    component = 0
    diffusivity = 'mu'
    face_variable = vel_bar_x
    prescribed_normal_flux = 0
    pressure_face_variable = pressure_bar
    pressure_variable = pressure
    variable = vel_x
    alpha = ${alpha}
  []
  [momentum_y_diffusion_neumann]
    type = NavierStokesStressIPHDGPrescribedTractionBC
    boundary = 'right'
    component = 1
    diffusivity = 'mu'
    face_variable = vel_bar_y
    prescribed_normal_flux = 0
    pressure_face_variable = pressure_bar
    pressure_variable = pressure
    variable = vel_y
    alpha = ${alpha}
  []
  [mass_neumann]
    type = MassContinuityIPHDGBC
    face_variable = pressure_bar
    variable = pressure
    boundary = 'right'
    face_velocity_functors = 'vel_bar_x vel_bar_y'
    interior_velocity_vars = 'vel_x vel_y'
  []
[]

[Materials]
  [const]
    type = ADGenericConstantMaterial
    prop_names = 'rho mu'
    prop_values = '${rho} ${mu}'
  []
  [vel]
    type = ADVectorFromComponentVariablesMaterial
    vector_prop_name = 'velocity'
    u = vel_x
    v = vel_y
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'
  nl_rel_tol = 1e-10
[]

[Postprocessors]
  [inlet_mass_flow]
    type = ConstantPostprocessor
    value = ${fparse -2*rho*Vin*width}
  []
#   [outlet_mass_flow]
#     type     = SideMassFluxIntegral
#     boundary = 'right'
#   []

  [Umax]
    type = ElementExtremeValue
    variable = vel_x
  []

#   [mass_balance_check]
#     type = FunctionValuePostprocessor
#     function = mass_balance
#   []
[]

# [Functions]
#   [mass_balance]
#     type = ParsedFunction
#     expression = '(mass_in + mass_out) / mass_out'
#     symbol_names = 'mass_in mass_out'
#     symbol_values = 'inlet_mass_flow outlet_mass_flow'
#   []
# []

[Outputs]
  exodus = true
  csv = true
  perf_graph = true
[]
