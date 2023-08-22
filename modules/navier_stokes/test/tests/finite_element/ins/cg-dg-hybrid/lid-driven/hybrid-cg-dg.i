mu = 1
rho = 1
l = 1
U = 1
n = 20

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${l}
    ymin = 0
    ymax = ${l}
    nx = ${n}
    ny = ${n}
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
    extra_matrix_tags = 'L'
  []
  [momentum_x_pressure]
    type = DGMomentumPressure
    variable = u
    pressure = pressure
    component = 0
  []
  [momentum_x_mass]
    type = Mass
    variable = u
    density = ${rho}
    matrix_tags = 'mass'
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
    extra_matrix_tags = 'L'
  []
  [momentum_y_pressure]
    type = DGMomentumPressure
    variable = v
    pressure = pressure
    component = 1
  []
  [momentum_y_mass]
    type = Mass
    variable = v
    density = ${rho}
    matrix_tags = 'mass'
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
    extra_matrix_tags = 'L'
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
    extra_matrix_tags = 'L'
  []
[]

[BCs]
  inactive = 'pressure_pin'
  [u_walls]
    type = DGFunctionDiffusionDirichletBC
    boundary = 'left bottom right'
    variable = u
    sigma = 6
    epsilon = -1
    function = '0'
    diff = 'mu'
    extra_matrix_tags = 'L'
  []
  [v_walls]
    type = DGFunctionDiffusionDirichletBC
    boundary = 'left bottom right top'
    variable = v
    sigma = 6
    epsilon = -1
    function = '0'
    diff = 'mu'
    extra_matrix_tags = 'L'
  []
  [u_top]
    type = DGFunctionDiffusionDirichletBC
    boundary = 'top'
    variable = u
    sigma = 6
    epsilon = -1
    function = '${U}'
    diff = 'mu'
    extra_matrix_tags = 'L'
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
    type = ProjectionAux
    variable = vel_x
    v = u
    execute_on = 'initial timestep_end'
  []
  [vel_y]
    type = ProjectionAux
    variable = vel_y
    v = v
    execute_on = 'initial timestep_end'
  []
  [p]
    type = ProjectionAux
    variable = p
    v = pressure
    execute_on = 'initial timestep_end'
  []
[]

[Problem]
  type = NavierStokesProblem
  mass_matrix = 'mass'
  L_matrix = 'L'
  extra_tag_matrices = 'mass L'
[]

[Preconditioning]
  active = FSP
  [FSP]
    type = FSP
    topsplit = 'up'
    [up]
      splitting = 'u p'
      splitting_type  = schur
      petsc_options_iname = '-pc_fieldsplit_schur_fact_type  -pc_fieldsplit_schur_precondition -ksp_gmres_restart -ksp_type -ksp_pc_side -ksp_rtol'
      petsc_options_value = 'full                            self                              300                fgmres    right        1e-4'
    []
    [u]
      vars = 'u v'
      # petsc_options = '-ksp_monitor'
      petsc_options_iname = '-pc_type -pc_hypre_type -ksp_type -ksp_rtol -ksp_gmres_restart -ksp_pc_side'
      petsc_options_value = 'hypre    boomeramg      gmres     1e-2      300                right'
    []
    [p]
      vars = 'pressure'
      petsc_options = '-pc_lsc_scale_diag -ksp_monitor'# -lsc_ksp_monitor'
      petsc_options_iname = '-ksp_type -ksp_gmres_restart -ksp_rtol -pc_type -ksp_pc_side -lsc_pc_type -lsc_ksp_type -lsc_ksp_pc_side -lsc_ksp_rtol'
      petsc_options_value = 'fgmres    300                1e-2      lsc      right        hypre        gmres         right            1e-1'
    []
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
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
