rho = 1000
mu = 1e-2
advected_interp_method='upwind'
velocity_interp_method='rc'
[Mesh]
  [mesh]
    type = CartesianMeshGenerator
    dim = 2
    dx = '1'
    ix = '10'
    dy = '1  1  5  1  1'
    iy = '10 10 50 10 10'
    subdomain_id = '3 2 1 2 4'
  []
[]
[GlobalParams]
  rhie_chow_user_object = 'rc'
[]
[UserObjects]
  [rc]
    type = PINSFVRhieChowInterpolator
    u = u_x
    v = u_y
    pressure = pressure
    porosity = porosity
    smoothing_layers = 100
    force_rc_correction = true
  []
[]
[Variables]
  [u_x]
    type = PINSFVSuperficialVelocityVariable
    initial_condition = 1e-6
  []
  [u_y]
    type = PINSFVSuperficialVelocityVariable
    initial_condition = 1e-6
  []
  [pressure]
    type = INSFVPressureVariable
  []
[]
[FVKernels]
  [mass]
    type = PINSFVMassAdvection
    variable = pressure
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    u = u_x
    v = u_y
    rho = ${rho}
    porosity = porosity
  []
  [u_time]
    type = INSFVMomentumTimeDerivative
    variable = u_x
    rho = ${rho}
    momentum_component = 'x'
  []
  [u_advection]
    type = PINSFVMomentumAdvection
    variable = u_x
    advected_quantity = 'rhou'
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    u = u_x
    v = u_y
    rho = ${rho}
    porosity = porosity
    momentum_component = 'x'
  []
  [u_pressure]
    type = PINSFVMomentumPressure
    variable = u_x
    pressure = pressure
    porosity = porosity
    momentum_component = 'x'
  []
  [v_time]
    type = INSFVMomentumTimeDerivative
    variable = u_y
    rho = ${rho}
    momentum_component = 'y'
  []
  [v_advection]
    type = PINSFVMomentumAdvection
    variable = u_y
    advected_quantity = 'rhov'
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    u = u_x
    v = u_y
    rho = ${rho}
    porosity = porosity
    momentum_component = 'y'
  []
  [v_pressure]
    type = PINSFVMomentumPressure
    variable = u_y
    pressure = pressure
    porosity = porosity
    momentum_component = 'y'
  []
  [u_drag]
    type = PINSFVMomentumFriction
    variable = u_x
    momentum_component = 'x'
    rho = ${rho}
    Darcy_name = darcy
  []
  [u_correct_oscillations]
    type = PINSFVOscillationsCorrection
    variable = u_x
    rho = ${rho}
    porosity = porosity
    Darcy_name = darcy
    consistent_scaling = 1e-2
    momentum_component = 'x'
  []
  [v_drag]
    type = PINSFVMomentumFriction
    variable = u_y
    momentum_component = 'y'
    rho = ${rho}
    Darcy_name = darcy
  []
  [v_correct_oscillations]
    type = PINSFVOscillationsCorrection
    variable = u_y
    rho = ${rho}
    porosity = porosity
    Darcy_name = darcy
    consistent_scaling = 1e-2
    momentum_component = 'y'
  []
  [u_viscosity]
    type = PINSFVMomentumDiffusion
    variable = u_x
    mu = ${mu}
    porosity = porosity
    momentum_component = 'x'
    superficial_velocity = 'velocity'
    #smooth_porosity = true
  []
  [v_viscosity]
    type = PINSFVMomentumDiffusion
    variable = u_y
    mu = ${mu}
    porosity = porosity
    momentum_component = 'y'
    superficial_velocity = 'velocity'
    #smooth_porosity = true
  []
[]
[AuxVariables]
  [eps_out]
    type = MooseVariableFVReal
  []
[]
[AuxKernels]
  [eps_out]
    type = ADFunctorElementalAux
    variable = eps_out
    functor = porosity
    execute_on = 'timestep_end'
  []
[]
[Functions]
  [ramp_v]
    type = PiecewiseLinear
    x = '0 10'
    y = '0 1'
  []
[]
[FVBCs]
  [inlet-u]
    type = INSFVInletVelocityBC
    boundary = 'bottom'
    variable = u_x
    function = 0
  []
  [inlet-v]
    type = INSFVInletVelocityBC
    boundary = 'bottom'
    variable = u_y
    function = ramp_v
  []
  [walls-u]
    type = INSFVNaturalFreeSlipBC
    boundary = 'left right'
    variable = u_x
    momentum_component = 'x'
  []
  [walls-v]
    type = INSFVNaturalFreeSlipBC
    boundary = 'left right'
    variable = u_y
    momentum_component = 'y'
  []
  [outlet_p]
    type = INSFVOutletPressureBC
    boundary = 'top'
    variable = pressure
    function = 0
  []
[]
[Materials]
  [ins_fv]
    type = INSFVMaterial
    u = 'u_x'
    v = 'u_y'
    pressure = 'pressure'
    rho = ${rho}
  []
  [Darcy_functor_bed]
    type = ADGenericVectorFunctorMaterial
    prop_names = 'darcy'
    prop_values = '200 200 200'
    block = '1'
  []
  [Darcy_functor]
    type = ADGenericVectorFunctorMaterial
    prop_names = 'darcy'
    prop_values = '50 50 50'
    block = '2 3 4'
  []
  [jump]
    type = ADPiecewiseConstantByBlockMaterial
    prop_name = 'porosity'
    subdomain_to_prop_value =  '
                               1 0.6
                               2 0.6
                               3 0.6
                               4 0.6
                               '
                               #'1 1 2 1 3 1 4 1 5 1'
  []
[]
[Executioner]
  type = Transient
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-2
    growth_factor = 1.1
    cutback_factor = 0.8
    iteration_window = 2
    optimal_iterations = 6
  []
  nl_max_its = 10
  end_time = 1e5
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'
  line_search = 'default'
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-8
[]
[Postprocessors]
  active = ''
  [inlet_p]
    type = SideAverageValue
    variable = 'pressure'
    boundary = 'inlet'
  []
  [outlet-u]
    type = SideIntegralVariablePostprocessor
    variable = u
    boundary = 'main_outlet'
  []
[]
[Outputs]
  exodus = true
  csv = true
[]
