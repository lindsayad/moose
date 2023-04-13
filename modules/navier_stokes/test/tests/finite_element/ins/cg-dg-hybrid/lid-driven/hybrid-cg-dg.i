mu = 1
rho = 1
l = 200
U = 1

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${l}
    ymin = 0
    ymax = ${l}
    nx = 20
    ny = 20
  []
  [corner_node]
    type = ExtraNodesetGenerator
    new_boundary = 'pinned_node'
    nodes = '0'
    input = gen
  []
[]

[Variables]
  [u]
    family = MONOMIAL
  []
  [v]
    family = MONOMIAL
  []
  [pressure][]
[]

[Kernels]
  [momentum_x_convection]
    type = ADConservativeAdvection
    variable = u
    velocity = 'velocity'
    advected_quantity = 'rhou'
  []
  [momentum_x_diffusion]
    type = MatDiffusion
    variable = u
    diffusivity = 'mu'
  []
  [momentum_x_pressure]
    type = DGMomentumPressure
    variable = u
    pressure = pressure
    component = 0
  []
  [momentum_y_convection]
    type = ADConservativeAdvection
    variable = v
    velocity = 'velocity'
    advected_quantity = 'rhov'
  []
  [momentum_y_diffusion]
    type = MatDiffusion
    variable = v
    diffusivity = 'mu'
  []
  [momentum_y_pressure]
    type = DGMomentumPressure
    variable = v
    pressure = pressure
    component = 1
  []
  [mass]
    type = CGMass
    variable = pressure
    velocity = velocity
  []
[]

[DGKernels]
  [momentum_x_convection]
    type = ADDGConvection
    variable = u
    velocity = 'velocity'
    advected_quantity = 'rhou'
  []
  [momentum_x_diffusion]
    type = DGDiffusion
    variable = u
    sigma = 6
    epsilon = -1
    diff = 'mu'
  []
  [momentum_y_convection]
    type = ADDGConvection
    variable = v
    velocity = 'velocity'
    advected_quantity = 'rhov'
  []
  [momentum_y_diffusion]
    type = DGDiffusion
    variable = v
    sigma = 6
    epsilon = -1
    diff = 'mu'
  []
[]

[BCs]
  [u_walls]
    type = DGFunctionDiffusionDirichletBC
    boundary = 'left bottom right'
    variable = u
    sigma = 6
    epsilon = -1
    function = '0'
    diff = 'mu'
  []
  [v_walls]
    type = DGFunctionDiffusionDirichletBC
    boundary = 'left bottom right top'
    variable = v
    sigma = 6
    epsilon = -1
    function = '0'
    diff = 'mu'
  []
  [u_top]
    type = DGFunctionDiffusionDirichletBC
    boundary = 'top'
    variable = u
    sigma = 6
    epsilon = -1
    function = '${U}'
    diff = 'mu'
  []
  [pressure_pin]
    type = DirichletBC
    variable = pressure
    boundary = 'pinned_node'
    value = 0
  []
[]

[Materials]
  [const]
    type = ADGenericConstantMaterial
    prop_names = 'rho'
    prop_values = '${rho}'
  []
  [const_reg]
    type = GenericConstantMaterial
    prop_names = 'mu'
    prop_values = '${mu}'
  []
  [vel]
    type = CGDGMaterial
    u = u
    v = v
    rho = 'rho'
  []
[]

[AuxVariables]
  [vel_x]
    family = MONOMIAL
    order = CONSTANT
  []
  [vel_y]
    family = MONOMIAL
    order = CONSTANT
  []
  [p]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [vel_x]
    type = SelfAux
    variable = vel_x
    v = u
  []
  [vel_y]
    type = SelfAux
    variable = vel_y
    v = v
  []
  [p]
    type = SelfAux
    variable = p
    v = pressure
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'
  nl_rel_tol = 1e-12
[]

[Outputs]
  exodus = true
[]

[Postprocessors]
  [Re]
    type = ParsedPostprocessor
    pp_names = ''
    function = '${rho} * ${U} * ${l} / ${mu}'
  []
[]
