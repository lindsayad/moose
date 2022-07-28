mu=1e-3
rho=1
velocity_interp_method = 'rc'
advected_interp_method = 'average'

[GlobalParams]
  rhie_chow_user_object = 'rc'
[]

[Mesh]
  [ccmg]
    type = ConcentricCircleMeshGenerator
    num_sectors = 6
    radii = '0.2546 0.3368'
    rings = '4 3 4'
    has_outer_square = on
    pitch = 1
    #portion = left_half
    preserve_volumes = off
    smoothing_max_it = 3
  []
  [rename_left]
    type = RenameBoundaryGenerator
    input = ccmg
    old_boundary = 'left'
    new_boundary = '101'
  []
  [left]
    type = CartesianMeshGenerator
    dim = 2
    dx = '5'
    dy = '1'
    ix = '100'
    iy = '16'
  []
  [move_it]
    type = TransformGenerator
    input = left
    transform = translate
    vector_value = '-5.5 -0.5 0'
  []
  [rename_right]
    type = RenameBoundaryGenerator
    input = move_it
    old_boundary = 'right'
    new_boundary = '102'
  []
  [stitch]
    type = StitchedMeshGenerator
    inputs = 'rename_left rename_right'
    stitch_boundaries_pairs = '101 102'
  []
  [in_between]
    type = SideSetsBetweenSubdomainsGenerator
    input = stitch
    primary_block = 2
    paired_block = 1
    new_boundary = 'no_circle'
  []
  [delete]
    type = BlockDeletionGenerator
    input = in_between
    block = '1'
  []
  [create_fused_top_sideset]
    input = delete
    type = ParsedGenerateSideset
    combinatorial_geometry = 'y > 0.49'
    normal = '0 1 0'
    new_sideset_name = 103
  []
  [top_left_block]
    type = GeneratedMeshGenerator
    xmin = -5.5
    xmax = -0.5
    ymin = 0.5
    ymax = ${fparse 0.5 + 2. / 16.}
    nx = 100
    ny = 2
    dim = 2
  []
  [rename_top_left_block]
    input = top_left_block
    type = RenameBlockGenerator
    old_block = '0'
    new_block = '100'
  []
  [rename_right_2]
    input = rename_top_left_block
    type = RenameBoundaryGenerator
    old_boundary = 'right'
    new_boundary = '104'
  []
  [top_right_block]
    type = GeneratedMeshGenerator
    xmin = -0.5
    xmax = 0.5
    ymin = 0.5
    ymax = ${fparse 0.5 + 2. / 16.}
    nx = 16
    ny = 2
    dim = 2
  []
  [rename_top_right_block]
    input = top_right_block
    type = RenameBlockGenerator
    old_block = '0'
    new_block = '101'
  []
  [rename_left_2]
    input = rename_top_right_block
    type = RenameBoundaryGenerator
    old_boundary = 'left'
    new_boundary = '105'
  []
  [stitch_2]
    inputs = 'rename_right_2 rename_left_2'
    type = StitchedMeshGenerator
    stitch_boundaries_pairs = '104 105'
  []
  [create_fused_bottom_sideset]
    input = stitch_2
    type = ParsedGenerateSideset
    combinatorial_geometry = 'y < 0.51'
    normal = '0 -1 0'
    new_sideset_name = 106
  []
  [stitch_3]
    inputs = 'create_fused_top_sideset create_fused_bottom_sideset'
    type = StitchedMeshGenerator
    stitch_boundaries_pairs = '103 106'
  []
  [no_slip_top]
    input = stitch_3
    type = ParsedGenerateSideset
    combinatorial_geometry = 'y > .615'
    normal = '0 1 0'
    new_sideset_name = 'no_slip_top'
  []
  [no_slip_bottom]
    input = no_slip_top
    type = ParsedGenerateSideset
    combinatorial_geometry = 'y < -0.49'
    normal = '0 -1 0'
    new_sideset_name = 'no_slip_bottom'
  []
  [inlet]
    input = no_slip_bottom
    type = ParsedGenerateSideset
    combinatorial_geometry = 'x > 0.49'
    normal = '1 0 0'
    new_sideset_name = 'inlet'
  []
  [outlet]
    input = inlet
    type = ParsedGenerateSideset
    combinatorial_geometry = 'x < -5.49'
    normal = '-1 0 0'
    new_sideset_name = 'outlet'
  []
  uniform_refine = 2
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = vel_x
    v = vel_y
    pressure = pressure
  []
[]

[Variables]
  [vel_x]
    type = INSFVVelocityVariable
    two_term_boundary_expansion = true
  []
  [vel_y]
    type = INSFVVelocityVariable
    two_term_boundary_expansion = true
  []
  [pressure]
    type = INSFVPressureVariable
    two_term_boundary_expansion = true
  []
[]

[FVKernels]
  # mass
  [mass]
    type = INSFVMassAdvection
    variable = pressure
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
  []
  [u_time]
    type = INSFVMomentumTimeDerivative
    variable = vel_x
    rho = ${rho}
    momentum_component = 'x'
  []
  [u_advection]
    type = INSFVMomentumAdvection
    variable = vel_x
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
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

  [v_time]
    type = INSFVMomentumTimeDerivative
    variable = vel_y
    rho = ${rho}
    momentum_component = 'y'
  []
  [v_advection]
    type = INSFVMomentumAdvection
    variable = vel_y
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
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
[]

[FVBCs]
  [no_slip_x]
    type = INSFVNoSlipWallBC
    variable = vel_x
    boundary = 'no_slip_top no_slip_bottom no_circle'
    function = 0
  []

  [no_slip_y]
    type = INSFVNoSlipWallBC
    variable = vel_y
    boundary = 'no_slip_top no_slip_bottom no_circle'
    function = 0
  []
  [inlet_x]
    type = INSFVInletVelocityBC
    variable = vel_x
    boundary = 'inlet'
    function = inlet_func
  []
  [inlet_y]
    type = INSFVInletVelocityBC
    variable = vel_y
    boundary = 'inlet'
    function = '1e-1 * exp(-t/10) * sin(t)'
  []
  [outlet_p]
    type = INSFVOutletPressureBC
    variable = pressure
    boundary = 'outlet'
    function = 0
  []
[]

[Functions]
  [inlet_func]
    type = ParsedFunction
    value = '-1'
  []
[]

[Materials]
  [const]
    type = GenericConstantMaterial
    prop_names = 'rho mu'
    prop_values = '${rho}  ${mu}'
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
    solve_type = 'NEWTON'
  []
[]

[Executioner]
  type = Transient
  dtmin = 1e-5
  petsc_options = '-snes_converged_reason -ksp_converged_reason'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  line_search = 'none'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-12
  nl_max_its = 10
  end_time = 50
  dtmax = 1
  scheme = 'bdf2'
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-5
    optimal_iterations = 6
    growth_factor = 1.5
  []
[]

[Outputs]
  exodus = true
  csv = true
  checkpoint = true
[]

[Postprocessors]
  [Re]
    type = ParsedPostprocessor
    function = 'rho * U * D / mu'
    constant_names = 'rho U D mu'
    constant_expressions = '${rho} 1 ${fparse 2 * .2546} ${mu}'
    pp_names = ''
  []
[]
