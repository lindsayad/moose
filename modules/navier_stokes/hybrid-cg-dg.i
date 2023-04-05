[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 100.0
    ymin = 0
    ymax = 100.0
    nx = 16
    ny = 16
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
  []
  [momentum_x_diffusion]
    type = Diffusion
    variable = u
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
  []
  [momentum_y_diffusion]
    type = Diffusion
    variable = v
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
  []
  [momentum_x_diffusion]
    type = DGDiffusion
    variable = u
    sigma = 6
    epsilon = -1
  []
  [momentum_y_convection]
    type = ADDGConvection
    variable = v
    velocity = 'velocity'
  []
  [momentum_y_diffusion]
    type = DGDiffusion
    variable = v
    sigma = 6
    epsilon = -1
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
  []
  [v_walls]
    type = DGFunctionDiffusionDirichletBC
    boundary = 'left bottom right top'
    variable = v
    sigma = 6
    epsilon = -1
    function = '0'
  []
  [u_top]
    type = DGFunctionDiffusionDirichletBC
    boundary = 'top'
    variable = u
    sigma = 6
    epsilon = -1
    function = '1'
  []
  [pressure_pin]
    type = DirichletBC
    variable = pressure
    boundary = 'pinned_node'
    value = 0
  []
[]

[Materials]
  [vel]
    type = CGDGMaterial
    u = u
    v = v
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
  [dof]
    type = DOFMap
    execute_on = 'initial'
  []
[]
