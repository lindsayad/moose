rho = 3279. # density [kg / m^3]  (@1000K)
mu = 0.005926 # viscosity [Pa s]
pump_force = '${fparse 1e-2 * 0.25*4.0e6}' # [N / m^3]

[Mesh]
  [gen]
    type = CartesianMeshGenerator
    dim = 2
    dx = '0.2 1.6 0.2'
    dy = '0.2 1.6 0.2'
    ix = '5 20 5'
    iy = '5 20 5'
    subdomain_id = '1 1 1
                    1 2 1
                    1 1 1'
  []
  [delete_internal_part]
    type = BlockDeletionGenerator
    input = gen
    block = '2'
    new_boundary = 'wall-internal'
  []
  [lump_bdries_to_wall]
    type = RenameBoundaryGenerator
    input = delete_internal_part
    old_boundary = 'bottom right top left'
    new_boundary = 'wall-external wall-external wall-external wall-external'
  []
  [pump_domain]
    type = ParsedSubdomainMeshGenerator
    input = lump_bdries_to_wall
    combinatorial_geometry = 'x > 0.6 & x < 0.8 & y > 1.0'
    block_id = '3'
  []
  [rename_blocks]
    type = RenameBlockGenerator
    input = pump_domain
    old_block = '1 3'
    new_block = 'pipe pump'
  []
  [side_pump]
    type = ParsedGenerateSideset
    input = rename_blocks
    included_subdomains = 'pump'
    included_neighbors = 'pipe'
    new_sideset_name = 'pump_side'
    normal = '1 0 0'
    combinatorial_geometry = 'x > 0.7'
  []
[]

[GlobalParams]
  velocity_interp_method = 'average'
  advected_interp_method = 'upwind'
  rhie_chow_user_object = 'rc'
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = vel_x
    v = vel_y
    pressure = pressure
    correct_volumetric_force = true
    volumetric_force_functors = 'pump_volume_force'
    volume_force_correction_method = 'pressure-consistent'
  []
[]

[Variables]
  [vel_x]
    type = INSFVVelocityVariable
    initial_condition = 1
  []
  [vel_y]
    type = INSFVVelocityVariable
    initial_condition = 1
  []
  [pressure]
    type = INSFVPressureVariable
  []
  # [lambda]
  #   family = SCALAR
  #   order = FIRST
  # []
[]

[AuxVariables]
  [U]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
  [Re]
    type = MooseVariableFVReal
  []
[]

[AuxKernels]
  [mag]
    type = VectorMagnitudeAux
    variable = U
    x = vel_x
    y = vel_y
  []
  [Re]
    type = ReynoldsNumberFunctorAux
    speed = U
    rho = ${rho}
    mu = ${mu}
    variable = Re
  []
[]

[FVKernels]
  [mass]
    type = INSFVMassAdvection
    variable = pressure
    rho = ${rho}
  []
  # [mean_zero_pressure]
  #   type = FVIntegralValueConstraint
  #   variable = pressure
  #   lambda = lambda
  #   phi0 = 0.0
  # []

  [u_time]
    type = INSFVMomentumTimeDerivative
    variable = vel_x
    rho = ${rho}
    momentum_component = 'x'
    extra_matrix_tags = 'mass'
  []
  [u_advection]
    type = INSFVMomentumAdvection
    variable = vel_x
    rho = ${rho}
    momentum_component = 'x'
  []
  [u_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_x
    mu = ${mu}
    momentum_component = 'x'
  []
  [u_pressure]
    type = INSFVMomentumPressure
    variable = vel_x
    momentum_component = 'x'
    pressure = pressure
  []
  [u_pump]
    type = INSFVBodyForce
    variable = vel_x
    functor = ${pump_force}
    block = 'pump'
    momentum_component = 'x'
  []
  # [u_mass]
  #   type = FVMass
  #   variable = vel_x
  #   matrix_tags = 'mass'
  #   density = ${rho}
  # []

  [v_time]
    type = INSFVMomentumTimeDerivative
    variable = vel_y
    rho = ${rho}
    momentum_component = 'y'
    extra_matrix_tags = 'mass'
  []
  [v_advection]
    type = INSFVMomentumAdvection
    variable = vel_y
    rho = ${rho}
    momentum_component = 'y'
  []
  [v_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_y
    mu = ${mu}
    momentum_component = 'y'
  []
  [v_pressure]
    type = INSFVMomentumPressure
    variable = vel_y
    momentum_component = 'y'
    pressure = pressure
  []
  # [v_mass]
  #   type = FVMass
  #   variable = vel_y
  #   matrix_tags = 'mass'
  #   density = ${rho}
  # []
[]

[FVBCs]
  [walls-u]
    type = INSFVNoSlipWallBC
    boundary = 'wall-internal wall-external'
    variable = vel_x
    function = '0'
  []
  [walls-v]
    type = INSFVNoSlipWallBC
    boundary = 'wall-internal wall-external'
    variable = vel_y
    function = '0'
  []
[]

[Functions]
  [pump_head]
    type = PiecewiseLinear
    x = '0.0 10.0'
    y = '1000.0 0.0'
  []
[]

[Materials]
  [pump_mat]
    type = NSFVPump
    rho = ${rho}
    speed = 'U'
    pressure_head_function = 'pump_head'
    rotation_speed = 120
    rotation_speed_rated = 100
    area_rated = 0.1
    volume_rated = 0.01
    flow_rate_rated = 1.0
    flow_rate = 'flow_rate'
    block = 'pump'
  []
[]

[Postprocessors]
  [flow_rate]
    type = Receiver
    default = 1.0
  []
  [flow_rate_to_pipe]
    type = VolumetricFlowRate
    advected_quantity = ${rho}
    boundary = 'pump_side'
    vel_x = 'vel_x'
    vel_y = 'vel_y'
  []
  [maximum_speed]
    type = ADElementExtremeFunctorValue
    functor = vel_x
    value_type = max
  []
  [maximum_cell_Re]
    type = ADElementExtremeFunctorValue
    functor = Re
    value_type = max
  []
[]

[Problem]
  extra_tag_matrices = 'mass L'
  type = NavierStokesProblem
  mass_matrix = 'mass'
  L_matrix = 'L'
  use_mass_matrix_for_scaling = false
  commute_lsc = false
  material_coverage_check = false
[]

[Preconditioning]
  [FSP]
    type = FSP
    topsplit = 'up'
    [up]
      splitting = 'u p'
      splitting_type  = schur
      petsc_options_iname = '-pc_fieldsplit_schur_fact_type  -pc_fieldsplit_schur_precondition -ksp_gmres_restart -ksp_rtol -ksp_type -ksp_atol -ksp_max_it'
      petsc_options_value = 'full                            self                             300                1e-5       fgmres    1e-9      10'
    []
      [u]
        vars = 'vel_x vel_y'
        petsc_options_iname = '-pc_type'
        petsc_options_value = 'lu'
      []
      [p]
        vars = 'pressure'
        petsc_options = '-ksp_converged_reason -pc_lsc_scale_diag'
        petsc_options_iname = '-ksp_type -ksp_gmres_restart -ksp_max_it -ksp_rtol -pc_type -ksp_pc_side -lsc_pc_type -lsc_pc_hypre_type -lsc_ksp_rtol -lsc_ksp_type -lsc_ksp_gmres_restart'
        petsc_options_value = 'fgmres    100                100         1e-2     lsc      right        hypre        boomeramg          1e-1          gmres         300'
      []
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-3
    optimal_iterations = 6
  []
  dtmax = 1
  steady_state_detection = true
  line_search = 'none'
  steady_state_tolerance = 1e-12
[]

[Outputs]
  exodus = true
  checkpoint = true
  [out]
    type = CSV
    execute_on = FINAL
    show = 'flow_rate_to_pipe maximum_speed'
  []
[]
