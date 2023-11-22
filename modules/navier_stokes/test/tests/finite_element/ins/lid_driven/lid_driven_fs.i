rho=1
mu=1
U=1
l=1
prefactor=${fparse 1/(l/2)^2}
n=2789

[GlobalParams]
  gravity = '0 0 0'
[]

[Mesh]
  [gen]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${l}
    ymin = 0
    ymax = ${l}
    nx = ${n}
    ny = ${n}
    elem_type = QUAD4
  []
  [./corner_node]
    type = ExtraNodesetGenerator
    new_boundary = 'pinned_node'
    nodes = '0'
    input = gen
  [../]
  second_order = true
  parallel_type = distributed
[]

[Variables]
  [./vel_x]
    order = SECOND
    family = LAGRANGE
  [../]
  [./vel_y]
    order = SECOND
    family = LAGRANGE
  [../]
  [./p]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./mass]
    type = INSMass
    variable = p
    u = vel_x
    v = vel_y
    pressure = p
  [../]

  [./x_momentum_space]
    type = INSMomentumLaplaceForm
    variable = vel_x
    u = vel_x
    v = vel_y
    pressure = p
    component = 0
    extra_matrix_tags = 'L'
  [../]
  [momentum_x_mass]
    type = MassKernel
    variable = vel_x
    density = ${rho}
    matrix_tags = 'mass'
  []
  [./y_momentum_space]
    type = INSMomentumLaplaceForm
    variable = vel_y
    u = vel_x
    v = vel_y
    pressure = p
    component = 1
    extra_matrix_tags = 'L'
  [../]
  [momentum_y_mass]
    type = MassKernel
    variable = vel_y
    density = ${rho}
    matrix_tags = 'mass'
  []
[]

[BCs]
  inactive = 'pressure_pin'
  [./x_no_slip]
    type = DirichletBC
    variable = vel_x
    boundary = 'bottom right left'
    value = 0.0
  [../]

  [./lid]
    type = FunctionDirichletBC
    variable = vel_x
    boundary = 'top'
    function = 'lid_function'
  [../]

  [./y_no_slip]
    type = DirichletBC
    variable = vel_y
    boundary = 'bottom right top left'
    value = 0.0
  [../]

  [./pressure_pin]
    type = DirichletBC
    variable = p
    boundary = 'pinned_node'
    value = 0
  [../]
[]

[Materials]
  [./const]
    type = GenericConstantMaterial
    block = 0
    prop_names = 'rho mu'
    prop_values = '${rho} ${mu}'
  [../]
[]

[Functions]
  [./lid_function]
    # We pick a function that is exactly represented in the velocity
    # space so that the Dirichlet conditions are the same regardless
    # of the mesh spacing.
    type = ParsedFunction
    expression = '${prefactor}*${U}*x*(${l}-x)'
  [../]
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
      vars = 'vel_x vel_y'
      petsc_options = '-ksp_monitor'
      petsc_options_iname = '-pc_type -pc_hypre_type -ksp_type -ksp_rtol -ksp_gmres_restart -ksp_pc_side'
      petsc_options_value = 'hypre    boomeramg      gmres     1e-2      300                right'
    []
    [p]
      vars = 'p'
      petsc_options = '-pc_lsc_scale_diag -ksp_monitor -lsc_ksp_monitor'
      petsc_options_iname = '-ksp_type -ksp_gmres_restart -ksp_rtol -pc_type -ksp_pc_side -lsc_pc_type -lsc_ksp_type -lsc_ksp_pc_side -lsc_ksp_rtol'
      petsc_options_value = 'fgmres    300                1e-2      lsc      right        hypre        gmres         right            1e-1'
    []
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
[]

[Outputs]
  exodus = false
  perf_graph = true
[]
