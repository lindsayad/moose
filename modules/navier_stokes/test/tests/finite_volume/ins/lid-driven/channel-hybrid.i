mu = 1.1
rho = 1.1

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 10
    ymin = -1
    ymax = 1
    nx = 100
    ny = 20
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

[Functions]
  [v_inlet]
    type = ParsedVectorFunction
    expression_x = '1'
  []
[]

[BCs]
  [u_walls]
    type = DGFunctionDiffusionDirichletBC
    boundary = 'bottom top'
    variable = u
    sigma = 6
    epsilon = -1
    function = '0'
    diff = 'mu'
  []
  [v_walls]
    type = DGFunctionDiffusionDirichletBC
    boundary = 'bottom top'
    variable = v
    sigma = 6
    epsilon = -1
    function = '0'
    diff = 'mu'
  []
  [u_in]
    type = ADDGConvectionBC
    boundary = 'left'
    variable = u
    velocity_function = v_inlet
    primal_dirichlet_value = 1
    rho = 'rho'
  []
  [v_in]
    type = ADDGConvectionBC
    boundary = 'left'
    variable = v
    velocity_function = v_inlet
    primal_dirichlet_value = 0
    rho = 'rho'
  []
  [p_in]
    type = ADDGConvectionBC
    boundary = 'left'
    variable = pressure
    velocity_function = v_inlet
    # Have to switch sign because we multiplied the mass equation weak
    # form by -1 in order to get a symmetric matrix
    primal_dirichlet_value = -1
  []
  [u_out]
    type = ADDGConvectionBC
    boundary = 'right'
    variable = u
    velocity_mat_prop = 'velocity'
    advected_quantity = 'rhou'
  []
  [v_out]
    type = ADDGConvectionBC
    boundary = 'right'
    variable = v
    velocity_mat_prop = 'velocity'
    advected_quantity = 'rhov'
  []
  [p_out]
    type = DirichletBC
    variable = pressure
    boundary = 'right'
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
