mu = 1.1
rho = 1.1
Vin = 1.0
gamma = 1e4
nu = '${fparse mu / rho}'

width = 1 # channel half-width
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
    ymin = '${fparse -width}'
    ymax = ${width}
    nx = 100
    ny = 20

    elem_type = QUAD4
  []
[]

[Problem]
  type = NavierStokesProblem
  extra_tag_matrices = 'mass'
  mass_matrix = 'mass'
  set_schur_pre = MASS
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
    family = LAGRANGE
    order = ${k}
  []
  [vel_bar_y]
    family = LAGRANGE
    order = ${k}
  []
  [pressure_bar]
    family = LAGRANGE
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

  [u_jump]
    type = MassFluxPenaltyIPHDG
    variable = vel_x
    face_variable = vel_bar_x
    face_velocity = face_velocity
    u = vel_x
    v = vel_y
    component = 0
    gamma = ${gamma}
  []
  [v_jump]
    type = MassFluxPenaltyIPHDG
    variable = vel_y
    face_variable = vel_bar_y
    face_velocity = face_velocity
    u = vel_x
    v = vel_y
    component = 1
    gamma = ${gamma}
  []
  [pb_mass]
    type = MassMatrixHDG
    variable = pressure_bar
    matrix_tags = 'mass'
    density = '${fparse -1/(gamma+nu)}'
  []
[]

[Functions]
  [u_inlet]
    type = ConstantFunction
    value = ${Vin}
  []
[]

[FunctorMaterials]
  [face_velocity]
    type = ADGenericVectorFunctorMaterial
    prop_names = face_velocity
    prop_values = 'vel_bar_x vel_bar_y 0'
  []
  [vel_inlet]
    type = GenericVectorFunctorMaterial
    prop_names = vel_inlet
    prop_values = 'u_inlet 0 0'
  []
  [vel_walls]
    type = GenericVectorFunctorMaterial
    prop_names = vel_walls
    prop_values = '0 0 0'
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

  [pb_mass]
    type = MassMatrixIntegratedBC
    variable = pressure_bar
    matrix_tags = 'mass'
    boundary = 'left right top bottom'
    density = '${fparse -1/(gamma+nu)}'
  []
  [u_jump_walls]
    type = MassFluxPenaltyBC
    variable = vel_x
    face_variable = vel_bar_x
    u = vel_x
    v = vel_y
    component = 0
    boundary = 'bottom top'
    gamma = ${gamma}
    face_velocity = vel_walls
    dirichlet_boundary = true
  []
  [v_jump_walls]
    type = MassFluxPenaltyBC
    variable = vel_y
    face_variable = vel_bar_y
    u = vel_x
    v = vel_y
    component = 1
    boundary = 'bottom top'
    gamma = ${gamma}
    face_velocity = vel_walls
    dirichlet_boundary = true
  []
  [u_jump_inlet]
    type = MassFluxPenaltyBC
    variable = vel_x
    face_variable = vel_bar_x
    u = vel_x
    v = vel_y
    component = 0
    boundary = 'left'
    gamma = ${gamma}
    face_velocity = vel_inlet
    dirichlet_boundary = true
  []
  [v_jump_inlet]
    type = MassFluxPenaltyBC
    variable = vel_y
    face_variable = vel_bar_y
    u = vel_x
    v = vel_y
    component = 1
    boundary = 'left'
    gamma = ${gamma}
    face_velocity = vel_inlet
    dirichlet_boundary = true
  []
  [u_jump_outlet]
    type = MassFluxPenaltyBC
    variable = vel_x
    face_variable = vel_bar_x
    u = vel_x
    v = vel_y
    component = 0
    boundary = 'right'
    gamma = ${gamma}
    face_velocity = face_velocity
    dirichlet_boundary = false
  []
  [v_jump_outlet]
    type = MassFluxPenaltyBC
    variable = vel_y
    face_variable = vel_bar_y
    u = vel_x
    v = vel_y
    component = 1
    boundary = 'right'
    gamma = ${gamma}
    face_velocity = face_velocity
    dirichlet_boundary = false
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

[Preconditioning]
  [FSP]
    type = SCFSP
    topsplit = 'up'
    [up]
      splitting = 'u p'
      splitting_type = schur
      petsc_options = '-ksp_converged_reason -ksp_monitor'
      petsc_options_iname = '-pc_fieldsplit_schur_fact_type  -pc_fieldsplit_schur_precondition -ksp_gmres_restart -ksp_type -ksp_pc_side -ksp_rtol -ksp_max_it -ksp_atol'
      petsc_options_value = 'lower                           self                              30                 fgmres    right        1e-4      30          1e-9'
    []
    [u]
      vars = 'vel_bar_x vel_bar_y'
      petsc_options = '-ksp_converged_reason'
      petsc_options_iname = '-pc_type -ksp_type -ksp_rtol -ksp_gmres_restart -ksp_pc_side -ksp_max_it -ksp_atol -ksp_norm_type   -pc_factor_mat_solver_type'
      petsc_options_value = 'lu       gmres     1e-2      300                right        300         1e-8      unpreconditioned mumps'
    []
    [p]
      vars = 'pressure_bar'
      petsc_options = '-ksp_ksp_converged_reason'
      petsc_options_iname = '-pc_type -ksp_pc_type -ksp_pc_jacobi_type -ksp_type -ksp_ksp_rtol -ksp_ksp_gmres_restart -ksp_ksp_pc_side -ksp_ksp_max_it -ksp_ksp_atol -ksp_ksp_norm_type'
      petsc_options_value = 'ksp      jacobi       diagonal            preonly   1e-2          300                    right            300             1e-8          unpreconditioned'
    []
  []
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-ksp_type'
  petsc_options_value = 'preonly'
  nl_rel_tol = 1e-10
  line_search = 'none'
[]

[Postprocessors]
  [inlet_mass_flow]
    type = ConstantPostprocessor
    value = '${fparse -2*rho*Vin*width}'
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
  print_linear_residuals = false
  exodus = true
  csv = true
  perf_graph = true
[]
