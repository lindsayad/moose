n=64
mu=1
rho=1

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

[AuxVariables]
  [vel_x]
    order = SECOND
  []
  [vel_y]
    order = SECOND
  []
[]

[AuxKernels]
  [vel_x]
    type = VectorVariableComponentAux
    variable = vel_x
    vector_variable = velocity
    component = 'x'
  []
  [vel_y]
    type = VectorVariableComponentAux
    variable = vel_y
    vector_variable = velocity
    component = 'y'
  []
[]

[Variables]
  [./velocity]
    order = SECOND
    family = LAGRANGE_VEC
  [../]

  [./p]
  [../]
[]

[Kernels]
  [./mass]
    type = INSADMass
    variable = p
  [../]
  [mass_mass]
    type = Mass
    variable = p
    matrix_tags = 'mass'
    density = ${rho}
  []

  # [./momentum_time]
  #   type = INSADMomentumTimeDerivative
  #   variable = velocity
  # [../]

  [./momentum_convection]
    type = INSADMomentumAdvection
    variable = velocity
  [../]

  [./momentum_viscous]
    type = INSADMomentumViscous
    variable = velocity
    extra_matrix_tags = 'L'
  [../]

  [./momentum_pressure]
    type = INSADMomentumPressure
    variable = velocity
    pressure = p
    integrate_p_by_parts = true
  [../]
  [velocity_mass]
    type = VectorMass
    variable = velocity
    matrix_tags = 'mass'
    density = ${rho}
  []
[]

[BCs]
  [./no_slip]
    type = VectorFunctionDirichletBC
    variable = velocity
    boundary = 'bottom right left'
    extra_matrix_tags = 'L'
  [../]

  [./lid]
    type = VectorFunctionDirichletBC
    variable = velocity
    boundary = 'top'
    function_x = 'lid_function'
    extra_matrix_tags = 'L'
  [../]
[]

[Materials]
  [./const]
    type = ADGenericConstantMaterial
    prop_names = 'rho mu'
    prop_values = '${rho} ${mu}'
  [../]
  [ins_mat]
    type = INSADMaterial
    velocity = velocity
    pressure = p
  []
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

[Problem]
  extra_tag_matrices = 'mass L'
  type = NavierStokesProblem
  mass_matrix = 'mass'
  L_matrix = 'L'
  use_pressure_mass_matrix = 'false'
  commute_lsc = 'true'
[]

[Preconditioning]
  active = 'FSP'
  [FSP]
    type = FSP
    topsplit = 'up'
    [up]
      splitting = 'u p'
      splitting_type  = schur
      petsc_options_iname = '-pc_fieldsplit_schur_fact_type  -pc_fieldsplit_schur_precondition -ksp_gmres_restart -ksp_rtol -ksp_type -ksp_atol'
      petsc_options_value = 'full                            self                              300                1e-5      fgmres     1e-8'
    []
      [u]
        vars = 'velocity'
        petsc_options_iname = '-pc_type -ksp_pc_side -ksp_type -ksp_rtol'
        petsc_options_value = 'lu       right        gmres     1e-5'
      []
      [p]
        vars = 'p'
        petsc_options = '-ksp_monitor'
        petsc_options_iname = '-ksp_type -ksp_gmres_restart -ksp_rtol -pc_type -ksp_pc_side -lsc_pc_type -lsc_mass_pc_type -lsc_mass_ksp_type -lsc_mass_ksp_pc_side'
        petsc_options_value = 'fgmres    300                1e-2      lsc      right        lu           lu                gmres              right'
      []
  []
[]

[Executioner]
  type = Steady
  line_search = 'none'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
[]

[Outputs]
  exodus = true
[]
