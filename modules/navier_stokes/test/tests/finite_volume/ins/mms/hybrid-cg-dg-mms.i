rho=1
mu=1

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = -1
    xmax = 1.0
    ymin = -1
    ymax = 1.0
    nx = 2
    ny = 2
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
    order = SECOND
  []
  [v]
    family = MONOMIAL
    order = SECOND
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
  [u_forcing]
    type = BodyForce
    variable = u
    function = forcing_u
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
  [v_forcing]
    type = BodyForce
    variable = v
    function = forcing_v
  []
  [mass]
    type = CGMass
    variable = pressure
    velocity = velocity
  []
  [p_forcing]
    type = BodyForce
    variable = pressure
    function = forcing_p
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
    boundary = 'left bottom right top'
    variable = u
    sigma = 6
    epsilon = -1
    function = exact_u
  []
  [v_walls]
    type = DGFunctionDiffusionDirichletBC
    boundary = 'left bottom right top'
    variable = v
    sigma = 6
    epsilon = -1
    function = exact_v
  []
  [pressure_pin]
    type = FunctionDirichletBC
    variable = pressure
    boundary = 'pinned_node'
    function = 'exact_p'
  []
[]

[Materials]
  [vel]
    type = CGDGMaterial
    u = u
    v = v
  []
[]

[Functions]
[exact_u]
  type = ParsedFunction
  value = 'sin(y)*cos((1/2)*x*pi)'
[]
[exact_rhou]
  type = ParsedFunction
  value = 'rho*sin(y)*cos((1/2)*x*pi)'
  vars = 'rho'
  vals = '${rho}'
[]
[forcing_u]
  type = ParsedFunction
  value = 'mu*sin(y)*cos((1/2)*x*pi) + (1/4)*pi^2*mu*sin(y)*cos((1/2)*x*pi) - 1/2*pi*rho*sin(x)*sin(y)*sin((1/2)*y*pi)*cos((1/2)*x*pi) + rho*sin(x)*cos(y)*cos((1/2)*x*pi)*cos((1/2)*y*pi) - pi*rho*sin(y)^2*sin((1/2)*x*pi)*cos((1/2)*x*pi) + sin(y)*cos(x)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[exact_v]
  type = ParsedFunction
  value = 'sin(x)*cos((1/2)*y*pi)'
[]
[exact_rhov]
  type = ParsedFunction
  value = 'rho*sin(x)*cos((1/2)*y*pi)'
  vars = 'rho'
  vals = '${rho}'
[]
[forcing_v]
  type = ParsedFunction
  value = 'mu*sin(x)*cos((1/2)*y*pi) + (1/4)*pi^2*mu*sin(x)*cos((1/2)*y*pi) - pi*rho*sin(x)^2*sin((1/2)*y*pi)*cos((1/2)*y*pi) - 1/2*pi*rho*sin(x)*sin(y)*sin((1/2)*x*pi)*cos((1/2)*y*pi) + rho*sin(y)*cos(x)*cos((1/2)*x*pi)*cos((1/2)*y*pi) + sin(x)*cos(y)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[exact_p]
  type = ParsedFunction
  value = 'sin(x)*sin(y)'
[]
[forcing_p]
  type = ParsedFunction
  value = '(1/2)*pi*rho*sin(x)*sin((1/2)*y*pi) + (1/2)*pi*rho*sin(y)*sin((1/2)*x*pi)'
  vars = 'rho'
  vals = '${rho}'
[]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       NONZERO               mumps'
  nl_rel_tol = 1e-12
[]

[Outputs]
  exodus = true
  csv = true
[]

[Postprocessors]
  [h]
    type = AverageElementSize
    outputs = 'console csv'
    execute_on = 'timestep_end'
  []
  [L2u]
    type = ElementL2Error
    variable = u
    function = exact_u
    outputs = 'console csv'
    execute_on = 'timestep_end'
  []
  [L2v]
    variable = v
    function = exact_v
    type = ElementL2Error
    outputs = 'console csv'
    execute_on = 'timestep_end'
  []
  [L2p]
    variable = pressure
    function = exact_p
    type = ElementL2Error
    outputs = 'console csv'
    execute_on = 'timestep_end'
  []
[]
