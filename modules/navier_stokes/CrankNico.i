# Mohammad Gomma, Mohamad Otoom, Abdullah Weiss
# Transient flow past a cylinder
# Testing
[GlobalParams]
  gravity = '0 0 0'
  integrate_p_by_parts = 'true'
[]

[Mesh]
  # Externally-generated mesh (GMSH)
  file = 111771_2ndOrder_Structured.msh
[]

[Variables]
  [vel_x]
    order = SECOND
    family = LAGRANGE
  []
  [vel_y]
    order = SECOND
    family = LAGRANGE
  []
  [p]
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  [x_momentum_time]
    type = INSMomentumTimeDerivative
    variable = vel_x
  []
  [y_momentum_time]
    type = INSMomentumTimeDerivative
    variable = vel_y
  []
  [mass]
    type = INSMass
    variable = p
    u = 'vel_x'
    v = 'vel_y'
    pressure = 'p'
  []
  [x_momentum_space]
    type = INSMomentumLaplaceForm
    variable = vel_x
    u = 'vel_x'
    v = 'vel_y'
    pressure = 'p'
    component = 0
  []
  [y_momentum_space]
    type = INSMomentumLaplaceForm
    variable = vel_y
    u = 'vel_x'
    v = 'vel_y'
    pressure = 'p'
    component = 1
  []
[]

[BCs]
  [x_no_slip]
    type = DirichletBC
    variable = vel_x
    boundary = 'top bottom wall'
    value = 0.0
  []
  [y_no_slip]
    type = DirichletBC
    variable = vel_y
    boundary = 'left top bottom wall'
    value = 0.0
  []
  [x_inlet]
    # function = 'inlet_func'
    type = DirichletBC
    variable = vel_x
    boundary = 'left'
    value = 0.2
  []
[]

[Materials]
  [const]
    type = GenericConstantMaterial
    block = '0'
    prop_names = 'rho mu'
    prop_values = '0.883  0.000883'
  []
[]

[Preconditioning]
  [SMP_PJFNK]
    # solve_type = PJFNK
    type = SMP
    full = true
    solve_type = NEWTON
  []
[]

[Executioner]
  type = Transient
  petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_levels -ksp_gmres_restart'
  petsc_options_value = 'asm      2               ilu          4                     300'
  line_search = none
  [./TimeIntegrator]
    type = CrankNicolson
  [../]
  dt = 0.02
  dtmin = 1E-6
  num_steps = 1500
  nl_rel_tol = 1e-8
  nl_max_its = 6
  l_tol = 1e-5
  l_max_its = 300
[]

[Outputs]
  exodus = true
  checkpoint = true
[]

# The following block generates csv file for each timestep for pressures around the cylinder wall (for drag and lift coefficients).
# Keep commented until the input file demonstrates good performance.
[VectorPostprocessors]
  [nodal_sample]
    # Pick off the wall pressure values.
    type = NodalValueSampler
    variable = 'p'
    boundary = 'wall'
    sort_by = x
  []
[]

[AuxVariables]
  [Co]
    family = MONOMIAL
    order = CONSTANT
  []
  [Re]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [Co]
    type = Courant
    execute_on = 'initial timestep_end'
    variable = Co
    u = vel_x
    v = vel_y
  []
  [Re]
    type = Reynolds
    execute_on = 'initial timestep_end'
    variable = Re
    u = vel_x
    v = vel_y
  []
[]
