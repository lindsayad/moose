[GlobalParams]
  displacements = 'disp_x disp_y'
  D_name = 1e3
  scaling = 1e0
  # use_displaced_mesh = true
[]

[Mesh]
  file = two-body-no-sep.e
[]

[Problem]
  kernel_coverage_check = false
[]

[Variables]
  [./disp_x]
    block = '1 2'
  [../]
  [./disp_y]
    block = '1 2'
  [../]
  [./lm]
    block = 3
  [../]
  [./tangent_lm]
    block = 3
  [../]
  [./vel_x]
    block = '1 2'
  [../]
  [./vel_y]
    block = '1 2'
  [../]
[]

[Kernels]
  [./disp_x]
    type = MatDiffusion
    variable = disp_x
  [../]
  [./disp_y]
    type = MatDiffusion
    variable = disp_y
  [../]
  [./accel_x]
    type = CoupledTimeDerivative
    variable = disp_x
    v = vel_x
    # block = 2
  [../]
  [./accel_y]
    type = CoupledTimeDerivative
    variable = disp_y
    v = vel_y
    # block = 2
  [../]
  [./coupled_time_velx]
    type = CoupledTimeDerivative
    variable = vel_x
    v = disp_x
    # block = 2
  [../]
  [./coupled_time_vely]
    type = CoupledTimeDerivative
    variable = vel_y
    v = disp_y
    # block = 2
  [../]
  [./source_velx]
    type = MatReaction
    variable = vel_x
    mob_name = 1
    # block = 2
  [../]
  [./source_vely]
    type = MatReaction
    variable = vel_y
    mob_name = 1
    # block = 2
  [../]
[]

# [Modules/TensorMechanics/Master]
#   [./all]
#     strain = SMALL
#     incremental = false
#     add_variables = true
#     generate_output = 'strain_xx strain_yy strain_zz' ## Not at all necessary, but nice
#     block = '1 2'
#   [../]
# []

# [Materials]
#   [./elasticity_tensor]
#     type = ComputeIsotropicElasticityTensor
#     youngs_modulus = 1e3
#     poissons_ratio = 0.3
#     block = '1 2'
#   [../]
#   [./small_stress]
#     type = ComputeLinearElasticStress
#     block = '1 2'
#   [../]
#   [./dummy]
#     type = GenericConstantMaterial
#     prop_names = 'dumb'
#     prop_values = '0'
#     block = 3
#   [../]
# []


[Constraints]
  [./lm]
    type = LMConstraint
    slave = 10
    master = 20
    variable = lm
    master_variable = disp_x
    disp_y = disp_y
  [../]
  [./tan_lm]
    type = TangentialLMConstraint
    lambda = 1
    slave = 10
    master = 20
    variable = tangent_lm
    contact_pressure = lm
    master_variable = vel_x
    vel_y = vel_y
    mu = 0.1
  [../]
[]

[BCs]
  [./botx]
    type = DirichletBC
    variable = disp_x
    boundary = 40
    value = 0.0
  [../]
  [./boty]
    type = DirichletBC
    variable = disp_y
    boundary = 40
    value = 0.0
  [../]
  [./topy]
    type = NeumannBC
    variable = disp_y
    boundary = 30
    value = -10e-4
  [../]
  [./leftx]
    type = NeumannBC
    variable = disp_x
    boundary = 50
    value = 5e-4
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 10
  dt = 10
  dtmin = 10
  solve_type = 'NEWTON'
  line_search = 'basic'
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_test_jacobian'# -snes_test_jacobian_view'
  petsc_options_iname = '-pc_type -snes_max_funcs'
  petsc_options_value = 'svd      100000'
  # nl_rel_tol = 1e-6

  l_max_its = 100
  nl_max_its = 1000
[]

[Outputs]
  exodus = true
  # checkpoint = true
  dofmap = true
[]

[Contact]
  [./leftright]
    master = 20
    slave = 10
    model = frictionless
    formulation = lagrange
    # penalty = 1e6
    system = constraint
    lm = lm
    tangent_lm = tangent_lm
    vel_x = vel_x
    vel_y = vel_y
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]
