n=64
mu=1

[GlobalParams]
  gravity = '0 0 0'
  preset = true
  supg = false
[]

[Problem]
  extra_tag_matrices = 'mass physics'
  previous_nl_solution_required = true
  type = NavierStokesProblem
  mass_matrix = 'mass'
  physics_matrix = 'physics'
  velocity_split_name = 'u'
  schur_fs_index = '1'
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 1.0
    ymin = 0
    ymax = 1.0
    nx = ${n}
    ny = ${n}
    elem_type = QUAD9
  []
  [./corner_node]
    type = ExtraNodesetGenerator
    new_boundary = 'pinned_node'
    nodes = '0'
    input = gen
  [../]
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
  # mass
  [./mass]
    type = INSMass
    variable = p
    u = vel_x
    v = vel_y
    pressure = p
    extra_matrix_tags = 'physics'
  [../]

  # x-momentum, space
  [./x_momentum_space]
    type = INSMomentumLaplaceForm
    variable = vel_x
    u = vel_x
    v = vel_y
    pressure = p
    component = 0
    extra_matrix_tags = 'physics'
  [../]
  [x_mass]
    type = Mass
    variable = vel_x
    matrix_tags = 'mass'
  []

  # y-momentum, space
  [./y_momentum_space]
    type = INSMomentumLaplaceForm
    variable = vel_y
    u = vel_x
    v = vel_y
    pressure = p
    component = 1
    extra_matrix_tags = 'physics'
  [../]
  [y_mass]
    type = Mass
    variable = vel_y
    matrix_tags = 'mass'
  []
[]

[BCs]
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

  # [./pressure_pin]
  #   type = DirichletBC
  #   variable = p
  #   boundary = 'pinned_node'
  #   value = 0
  # [../]
[]

[Materials]
  [./const]
    type = GenericConstantMaterial
    block = 0
    prop_names = 'rho mu'
    prop_values = '1  ${mu}'
  [../]
[]

[Functions]
  [./lid_function]
    # We pick a function that is exactly represented in the velocity
    # space so that the Dirichlet conditions are the same regardless
    # of the mesh spacing.
    type = ParsedFunction
    expression = '4*x*(1-x)'
  [../]
[]

[Preconditioning]
  active = 'FSP'
  [FSP]
    type = FSP
    topsplit = 'by_diri_others'
    [by_diri_others]
      splitting = 'diri others'
      splitting_type  = multiplicative
      petsc_options_iname = '-ksp_rtol -ksp_type -ksp_pc_side'
      petsc_options_value = '1e-5      fgmres    right'
    []
      [diri]
        sides = 'left right top bottom'
        vars = 'vel_x vel_y'
        petsc_options_iname = '-pc_type -ksp_pc_side -ksp_type -ksp_rtol -pc_factor_mat_solver_type'
        petsc_options_value = 'lu       right        gmres     1e-5      mumps'
      []
      [others]
        petsc_options_iname = '-pc_fieldsplit_schur_fact_type  -pc_fieldsplit_schur_precondition -ksp_gmres_restart -ksp_rtol -ksp_type -ksp_atol -ksp_rtol -ksp_pc_side'
        petsc_options_value = 'full                            self                             300                1e-5      fgmres     1e-9      1e-5      right'
        splitting = 'u p'
        splitting_type = schur
        unside_by_var_boundary_name = 'left top right bottom left top right bottom'
        unside_by_var_var_name = 'vel_x vel_x vel_x vel_x vel_y vel_y vel_y vel_y'
        petsc_options = '-ksp_monitor'
      []
        [u]
          unside_by_var_boundary_name = 'left top right bottom left top right bottom'
          unside_by_var_var_name = 'vel_x vel_x vel_x vel_x vel_y vel_y vel_y vel_y'
          vars = 'vel_x vel_y'
          petsc_options_iname = '-pc_type -ksp_pc_side -ksp_type -ksp_rtol -pc_factor_mat_solver_type'
          petsc_options_value = 'lu       right        gmres     1e-5      mumps'
        []
        [p]
          vars = 'p'
          petsc_options = '-ksp_monitor -pc_lsc_scale_diag'
          petsc_options_iname = '-ksp_type -ksp_gmres_restart -ksp_rtol -pc_type -ksp_pc_side -pc_type  -lsc_pc_type -lsc_pc_factor_mat_solver_type -lsc_ksp_type -lsc_ksp_rtol -lsc_ksp_pc_side -lsc_ksp_gmres_restart'
          petsc_options_value = 'fgmres     300                1e-5     lsc      right        lsc       lu           mumps                          gmres         1e-5          right            300'
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
  solve_type = NEWTON
  type = Steady
  petsc_options_iname = '-snes_max_it'
  petsc_options_value = '100'
  line_search = 'none'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
[]

[Outputs]
  exodus = true
[]
