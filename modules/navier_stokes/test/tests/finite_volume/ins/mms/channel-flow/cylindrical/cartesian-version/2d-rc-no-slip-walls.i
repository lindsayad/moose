mu=1.1
rho=1.1

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = -1
    xmax = 1
    ymin = -1
    ymax = 1
    nx = 2
    ny = 2
  []
[]

[Problem]
  fv_bcs_integrity_check = false
[]

[GlobalParams]
  rhie_chow_user_object = 'rc'
  two_term_boundary_expansion = true
  advected_interp_method = 'average'
  velocity_interp_method = 'rc'
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = u
    v = v
    standard_body_forces = true
  []
[]

[Variables]
  [u]
    type = INSFVVelocityVariable
    # initial_condition = 1e-15
  []
  [v]
    type = INSFVVelocityVariable
    # initial_condition = 1e-15
  []
  [pressure]
    type = INSFVPressureVariable
  []
  # [lambda]
  #   family = SCALAR
  #   order = FIRST
  # []
[]

[ICs]
  [u]
    type = FunctionIC
    function = 'exact_u'
    variable = u
  []
  [v]
    type = FunctionIC
    function = 'exact_v'
    variable = v
  []
  [pressure]
    type = FunctionIC
    function = 'exact_p'
    variable = pressure
  []
[]

[FVKernels]
  [mass]
    type = INSFVMassAdvection
    variable = pressure
    vel = 'velocity'
    pressure = pressure
    u = u
    v = v
    rho = ${rho}
    # boundaries_to_force = 'top'
  []
  [mass_forcing]
    type = FVBodyForce
    variable = pressure
    function = forcing_p
  []
  # [mean_zero_pressure]
  #   type = FVScalarLagrangeMultiplier
  #   variable = pressure
  #   lambda = lambda
  # []

  [u_advection]
    type = INSFVMomentumAdvection
    variable = u
    advected_quantity = 'rhou'
    vel = 'velocity'
    pressure = pressure
    u = u
    v = v
    rho = ${rho}
    momentum_component = 'x'
    # boundaries_to_force = 'top'
  []
  [u_viscosity]
    type = INSFVMomentumDiffusion
    variable = u
    mu = ${mu}
    momentum_component = 'x'
  []
  [u_pressure]
    type = INSFVMomentumPressure
    variable = u
    momentum_component = 'x'
    pressure = pressure
  []
  [u_forcing]
    type = INSFVBodyForce
    variable = u
    function = forcing_u
    momentum_component = 'x'
  []

  [v_advection]
    type = INSFVMomentumAdvection
    variable = v
    advected_quantity = 'rhov'
    vel = 'velocity'
    pressure = pressure
    u = u
    v = v
    rho = ${rho}
    momentum_component = 'y'
    # boundaries_to_force = 'top'
  []
  [v_viscosity]
    type = INSFVMomentumDiffusion
    variable = v
    mu = ${mu}
    momentum_component = 'y'
  []
  [v_pressure]
    type = INSFVMomentumPressure
    variable = v
    momentum_component = 'y'
    pressure = pressure
  []
  [v_forcing]
    type = INSFVBodyForce
    variable = v
    function = forcing_v
    momentum_component = 'y'
  []
[]

[FVBCs]
  [top_u_diffusive_flux]
    type = INSFVMomentumFunctionFluxBC
    variable = u
    boundary = 'top'
    function = 'diffusive_flux_u_top'
    momentum_component = 'x'
  []
  [top_v_diffusive_flux]
    type = INSFVMomentumFunctionFluxBC
    variable = v
    boundary = 'top'
    function = 'diffusive_flux_v_top'
    momentum_component = 'y'
  []
  [left_pressure_flux]
    type = FVFunctionNeumannBC
    variable = pressure
    boundary = 'left'
    function = 'flux_p_left'
  []
  [right_pressure_flux]
    type = FVFunctionNeumannBC
    variable = pressure
    boundary = 'right'
    function = 'flux_p_right'
  []
# Dirichlet conditions for velocity
  [u_walls]
    type = INSFVNoSlipWallBC
    variable = u
    boundary = 'left right'
    function = 'exact_u'
  []
  [v_walls]
    type = INSFVNoSlipWallBC
    variable = v
    boundary = 'left right'
    function = 'exact_v'
  []
  # [p]
  #   type = FVFunctionNeumannBC
  #   variable = pressure
  #   boundary = 'top'
  #   function = 'flux_p_top'
  # []
  [p]
    type = FVFunctionDirichletBC
    variable = pressure
    function = 'exact_p'
    boundary = 'top'
  []
  # [p]
  #   type = INSFVOutletPressureBC
  #   variable = pressure
  #   function = 'exact_p'
  #   boundary = 'top'
  # []
  # [p]
  #   type = INSFVMassAdvectionOutflowBC
  #   variable = pressure
  #   boundary = 'top'
  #   rho = ${rho}
  #   u = u
  #   v = v
  # []
  [u_advection_top]
    type = INSFVMomentumAdvectionOutflowBC
    variable = u
    boundary = 'top'
    u = u
    v = v
    rho = ${rho}
    momentum_component = 'x'
    advected_quantity = 'rhou'
    vel = 'velocity'
  []
  [v_advection_top]
    type = INSFVMomentumAdvectionOutflowBC
    variable = v
    boundary = 'top'
    u = u
    v = v
    rho = ${rho}
    momentum_component = 'y'
    advected_quantity = 'rhov'
    vel = 'velocity'
  []
  # # These are dirichlet boundary conditions in the background that trigger
  # # execution of the kernels on the boundary
  [inlet_u]
    type = INSFVInletVelocityBC
    variable = u
    function = 'exact_u'
    boundary = 'bottom'
  []
  [inlet_v]
    type = INSFVInletVelocityBC
    variable = v
    function = 'exact_v'
    boundary = 'bottom'
  []
[]

[Materials]
  [ins_fv]
    type = INSFVMaterial
    u = 'u'
    v = 'v'
    pressure = 'pressure'
    rho = ${rho}
  []
[]

[Functions]
[exact_u]
  type = ParsedFunction
  value = 'sin(1.1*x)*cos(1.2*y)'
[]
[forcing_u]
  type = ParsedFunction
  value = '2.65*mu*sin(1.1*x)*cos(1.2*y) - 1.2*rho*sin(1.1*x)*sin(1.2*y)*cos(1.3*x)*cos(1.4*y) - 1.4*rho*sin(1.1*x)*sin(1.4*y)*cos(1.3*x)*cos(1.2*y) + 2.2*rho*sin(1.1*x)*cos(1.1*x)*cos(1.2*y)^2 + 1.5*cos(1.5*x)*cos(1.6*y)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[advective_flux_u_left]
  type = ParsedFunction
  value = '-rho*sin(1.1*x)^2*cos(1.2*y)^2'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[advective_flux_u_right]
  type = ParsedFunction
  value = 'rho*sin(1.1*x)^2*cos(1.2*y)^2'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[advective_flux_u_top]
  type = ParsedFunction
  value = 'rho*sin(1.1*x)*cos(1.3*x)*cos(1.2*y)*cos(1.4*y)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[advective_flux_u_bottom]
  type = ParsedFunction
  value = '-rho*sin(1.1*x)*cos(1.3*x)*cos(1.2*y)*cos(1.4*y)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[diffusive_flux_u_left]
  type = ParsedFunction
  value = '1.1*mu*cos(1.1*x)*cos(1.2*y)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[diffusive_flux_u_right]
  type = ParsedFunction
  value = '-1.1*mu*cos(1.1*x)*cos(1.2*y)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[diffusive_flux_u_top]
  type = ParsedFunction
  value = '1.2*mu*sin(1.1*x)*sin(1.2*y)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[diffusive_flux_u_bottom]
  type = ParsedFunction
  value = '-1.2*mu*sin(1.1*x)*sin(1.2*y)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[exact_v]
  type = ParsedFunction
  value = 'cos(1.3*x)*cos(1.4*y)'
[]
[forcing_v]
  type = ParsedFunction
  value = '3.65*mu*cos(1.3*x)*cos(1.4*y) - 1.3*rho*sin(1.1*x)*sin(1.3*x)*cos(1.2*y)*cos(1.4*y) - 2.8*rho*sin(1.4*y)*cos(1.3*x)^2*cos(1.4*y) + 1.1*rho*cos(1.1*x)*cos(1.3*x)*cos(1.2*y)*cos(1.4*y) - 1.6*sin(1.5*x)*sin(1.6*y)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[advective_flux_v_left]
  type = ParsedFunction
  value = '-rho*sin(1.1*x)*cos(1.3*x)*cos(1.2*y)*cos(1.4*y)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[advective_flux_v_right]
  type = ParsedFunction
  value = 'rho*sin(1.1*x)*cos(1.3*x)*cos(1.2*y)*cos(1.4*y)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[advective_flux_v_top]
  type = ParsedFunction
  value = 'rho*cos(1.3*x)^2*cos(1.4*y)^2'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[advective_flux_v_bottom]
  type = ParsedFunction
  value = '-rho*cos(1.3*x)^2*cos(1.4*y)^2'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[diffusive_flux_v_left]
  type = ParsedFunction
  value = '-1.3*mu*sin(1.3*x)*cos(1.4*y)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[diffusive_flux_v_right]
  type = ParsedFunction
  value = '1.3*mu*sin(1.3*x)*cos(1.4*y)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[diffusive_flux_v_top]
  type = ParsedFunction
  value = '1.4*mu*sin(1.4*y)*cos(1.3*x)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[diffusive_flux_v_bottom]
  type = ParsedFunction
  value = '-1.4*mu*sin(1.4*y)*cos(1.3*x)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[exact_p]
  type = ParsedFunction
  value = 'sin(1.5*x)*cos(1.6*y)'
[]
[forcing_p]
  type = ParsedFunction
  value = '-1.4*rho*sin(1.4*y)*cos(1.3*x) + 1.1*rho*cos(1.1*x)*cos(1.2*y)'
  vars = 'rho'
  vals = '${rho}'
[]
[flux_p_left]
  type = ParsedFunction
  value = '-rho*sin(1.1*x)*cos(1.2*y)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[flux_p_right]
  type = ParsedFunction
  value = 'rho*sin(1.1*x)*cos(1.2*y)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[flux_p_top]
  type = ParsedFunction
  value = 'rho*cos(1.3*x)*cos(1.4*y)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[flux_p_bottom]
  type = ParsedFunction
  value = '-rho*cos(1.3*x)*cos(1.4*y)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       NONZERO               superlu_dist'
  line_search = 'none'
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-12
[]

[Outputs]
  exodus = true
  csv = true
  [dof]
    type = DOFMap
    execute_on = 'initial'
  []
[]

[Postprocessors]
  [h]
    type = AverageElementSize
    outputs = 'console csv'
    execute_on = 'timestep_end'
  []
  [./L2u]
    type = ElementL2Error
    variable = u
    function = exact_u
    outputs = 'console csv'
    execute_on = 'timestep_end'
  [../]
  [./L2v]
    type = ElementL2Error
    variable = v
    function = exact_v
    outputs = 'console csv'
    execute_on = 'timestep_end'
  [../]
  [./L2p]
    variable = pressure
    function = exact_p
    type = ElementL2Error
    outputs = 'console csv'
    execute_on = 'timestep_end'
  [../]
[]
