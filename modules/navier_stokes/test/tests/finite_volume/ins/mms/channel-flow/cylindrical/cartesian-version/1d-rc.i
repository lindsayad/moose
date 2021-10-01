mu=1
rho=1
advected_interp_method='average'
velocity_interp_method='rc'

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = -1
    xmax = 1
    nx = 2
  []
[]

[Problem]
  fv_bcs_integrity_check = false
[]

[GlobalParams]
  rhie_chow_user_object = 'rc'
  two_term_boundary_expansion = true
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = u
    standard_body_forces = true
  []
[]

[Variables]
  [u]
    type = INSFVVelocityVariable
    initial_condition = 1e-15
  []
  [pressure]
    type = INSFVPressureVariable
  []
  [lambda]
    family = SCALAR
    order = FIRST
  []
[]

[FVKernels]
  [mass]
    type = INSFVMassAdvection
    variable = pressure
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    vel = 'velocity'
    pressure = pressure
    u = u
    rho = ${rho}
  []
  [mass_forcing]
    type = FVBodyForce
    variable = pressure
    function = forcing_p
  []
  [mean_zero_pressure]
    type = FVScalarLagrangeMultiplier
    variable = pressure
    lambda = lambda
  []

  [u_advection]
    type = INSFVMomentumAdvection
    variable = u
    advected_quantity = 'rhou'
    vel = 'velocity'
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    pressure = pressure
    u = u
    rho = ${rho}
    momentum_component = 'x'
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
[]

[FVBCs]
  # We envision an infinite channel flow problem with symmetry boundary on the
  # left and something like a free slip BC on the right; critically these are
  # both flux bcs.  Flux BC in production means to mimic it in MMS we can add
  # arbitrarily many FVFunctionNeumannBCs and/or INSFVMomentumFunctionFluxBCs to
  # get the boundary flux correct

  # This corresponds exactly to a diffusive flux for the velocity component that
  # is perpendicular to the boundary
  [symmetry]
    type = INSFVSymmetryVelocityBC
    boundary = 'left'
    variable = u
    u = u
    momentum_component = 'x'
    mu = ${mu}
  []

  # And since we are doing MMS we need to toss in the advective flux at the left
  # boundary
  [u_advection_left]
    type = INSFVMomentumFunctionFluxBC
    variable = u
    boundary = 'left'
    function = 'flux_u_left'
    momentum_component = 'x'
  []

  # Free slip flux bc right which is 0 flux so we add MMS bcs for both advective
  # and diffusive flux contributions
  [u_advection_right]
    type = INSFVMomentumFunctionFluxBC
    variable = u
    boundary = 'right'
    function = 'flux_u_right'
    momentum_component = 'x'
  []
  [u_diffusion_right]
    type = INSFVMomentumFunctionFluxBC
    variable = u
    boundary = 'right'
    function = 'flux_u_diffusion_right'
    momentum_component = 'x'
  []

  # Prescribe fluxes. In our idealized case we would be applying natural
  # boundary conditions (which are really just 0 flux bcs) to the mass equation
  # on both left and right so that means we add MMS flux BCs for the mass
  # advective flux
  [pressure_left]
    type = FVFunctionNeumannBC
    variable = pressure
    boundary = 'left'
    function = 'flux_p_left'
  []
  [pressure_right]
    type = FVFunctionNeumannBC
    variable = pressure
    boundary = 'right'
    function = 'flux_p_right'
  []

  # Finally a symmetry bc means we know that the perpendicular velocity is zero,
  # e.g. we know a Dirichlet condition for it, so we correspondingly apply an
  # MMS dirichlet BC here
  [diri_u]
    type = FVFunctionDirichletBC
    variable = u
    function = 'exact_u'
    boundary = 'left'
  []
[]

[Materials]
  [ins_fv]
    type = INSFVMaterial
    u = 'u'
    pressure = 'pressure'
    rho = ${rho}
  []
[]

[Functions]
[exact_u]
  type = ParsedFunction
  value = 'cos(x)'
[]
[forcing_u]
  type = ParsedFunction
  value = 'mu*cos(x) - 2*rho*sin(x)*cos(x) + cos(x)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[flux_u_left]
  type = ParsedFunction
  value = '-rho*cos(x)^2'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[flux_u_right]
  type = ParsedFunction
  value = 'rho*cos(x)^2'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[flux_u_diffusion_left]
  type = ParsedFunction
  value = '-mu*sin(x)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[flux_u_diffusion_right]
  type = ParsedFunction
  value = 'mu*sin(x)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[exact_p]
  type = ParsedFunction
  value = 'sin(x)'
[]
[forcing_p]
  type = ParsedFunction
  value = '-rho*sin(x)'
  vars = 'rho'
  vals = '${rho}'
[]
[flux_p_left]
  type = ParsedFunction
  value = '-rho*cos(x)'
  vars = 'rho'
  vals = '${rho}'
[]
[flux_p_right]
  type = ParsedFunction
  value = 'rho*cos(x)'
  vars = 'rho'
  vals = '${rho}'
[]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      200                lu           NONZERO'
  line_search = 'none'
  nl_rel_tol = 1e-12
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
  [./L2p]
    variable = pressure
    function = exact_p
    type = ElementL2Error
    outputs = 'console csv'
    execute_on = 'timestep_end'
  [../]
[]
