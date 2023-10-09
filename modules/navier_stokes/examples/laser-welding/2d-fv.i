period=.35e-4 # s
endtime=${fparse 3 * period} # s
timestep=${fparse period / 100} # s
surfacetemp=2700 # K
bottomtemp=2700 # K
sb=5.67e-8 # W/(m^2 K^4)
advected_interp_method='upwind'
velocity_interp_method='rc'
rho='rho'
mu='mu'
cp='cp'

[GlobalParams]
  rhie_chow_user_object = 'rc'
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = -.7e-3 # m
  xmax = 0.7e-3 # m
  ymin = -.35e-3 # m
  ymax = 0
  nx = 75
  ny = 20
  displacements = 'disp_x disp_y'
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = vel_x
    v = vel_y
    pressure = pressure
    use_displaced_mesh = true
    disp_x = disp_x
    disp_y = disp_y
  []
[]

[AuxVariables]
  [mu_out]
    type = MooseVariableFVReal
  []
[]

[AuxKernels]
  [mu_out]
    type = FunctorElementalAux
    functor = mu
    variable = mu_out
    execute_on = timestep_end
  []
[]

[Variables]
  [vel_x]
    type = INSFVVelocityVariable
  []
  [vel_y]
    type = INSFVVelocityVariable
  []
  [T]
    type = INSFVEnergyVariable
  []
  [pressure]
    type = INSFVPressureVariable
  []
  [disp_x]
  []
  [disp_y]
  []
[]

[ICs]
  [T]
    type = FunctionIC
    variable = T
    function = '2700 + ((${surfacetemp} - ${bottomtemp}) / .35e-3) * y'
  []
[]

[Kernels]
  [disp_x]
    type = MatDiffusion
    variable = disp_x
    diffusivity = 1e6
  []
  [disp_y]
    type = MatDiffusion
    variable = disp_y
    diffusivity = 1e6
  []
[]

[FVKernels]
  # pressure equation
  [mass]
    type = INSFVMassAdvection
    variable = pressure
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
    use_displaced_mesh = true
    boundaries_to_force = top
  []

  # momentum equations
  # u equation
  [u_time]
    type = INSFVMomentumTimeDerivative
    variable = vel_x
    rho = ${rho}
    momentum_component = 'x'
    use_displaced_mesh = true
  []
  [u_advection]
    type = INSFVMomentumAdvection
    variable = vel_x
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
    momentum_component = 'x'
    use_displaced_mesh = true
  []
  [u_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_x
    mu = ${mu}
    momentum_component = 'x'
    use_displaced_mesh = true
  []
  [u_pressure]
    type = INSFVMomentumPressureFlux
    variable = vel_x
    momentum_component = 'x'
    pressure = pressure
    use_displaced_mesh = true
  []
  # v equation
  [v_time]
    type = INSFVMomentumTimeDerivative
    variable = vel_y
    rho = ${rho}
    momentum_component = 'y'
    use_displaced_mesh = true
  []
  [v_advection]
    type = INSFVMomentumAdvection
    variable = vel_y
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
    momentum_component = 'y'
    use_displaced_mesh = true
  []
  [v_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_y
    mu = ${mu}
    momentum_component = 'y'
    use_displaced_mesh = true
  []
  [v_pressure]
    type = INSFVMomentumPressureFlux
    variable = vel_y
    momentum_component = 'y'
    pressure = pressure
    use_displaced_mesh = true
  []
  # energy equation
  [temperature_time]
    type = INSFVEnergyTimeDerivative
    variable = T
    rho = ${rho}
    cp = ${cp}
    use_displaced_mesh = true
  []
  [temperature_advection]
    type = INSFVEnergyAdvection
    variable = T
    use_displaced_mesh = true
  []
  [temperature_conduction]
    type = FVDiffusion
    coeff = 'k'
    variable = T
    use_displaced_mesh = true
  []
[]

[FVBCs]
  # momentum boundary conditions
  [no_slip_x]
    type = INSFVNoSlipWallBC
    variable = vel_x
    boundary = 'bottom right left'
    function = 0
  []
  [no_slip_y]
    type = INSFVNoSlipWallBC
    variable = vel_y
    boundary = 'bottom right left'
    function = 0
  []
  [vapor_recoil_x]
    type = INSFVVaporRecoilPressureMomentumFluxBC
    variable = vel_x
    boundary = 'top'
    momentum_component = 'x'
    rc_pressure = rc_pressure
    use_displaced_mesh = true
  []
  [vapor_recoil_y]
    type = INSFVVaporRecoilPressureMomentumFluxBC
    variable = vel_y
    boundary = 'top'
    momentum_component = 'y'
    rc_pressure = rc_pressure
    use_displaced_mesh = true
  []
  # energy boundary conditions
  [T_cold]
    type = FVDirichletBC
    variable = T
    boundary = 'bottom'
    value = '${bottomtemp}'
  []
  [radiation_flux]
    type = FVFunctorRadiativeBC
    variable = T
    boundary = 'top'
    emissivity = '1'
    Tinfinity = 300
    stefan_boltzmann_constant = ${sb}
    use_displaced_mesh = true
  []
  [weld_flux]
    type = FVGaussianEnergyFluxBC
    variable = T
    boundary = 'top'
    P0 = 159.96989792079225
    R = 1.25e-4
    x_beam_coord = '2e-4 * sin(t * 2 * pi / ${period})'
    y_beam_coord = 0
    z_beam_coord = 0
    use_displaced_mesh = true
  []
[]

[BCs]
  # displacement boundary conditions
  [x_no_disp]
    type = DirichletBC
    variable = disp_x
    boundary = 'bottom'
    value = 0
  []
  [y_no_disp]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  []
  [displace_x_top]
    type = INSADDisplaceBoundaryBC
    boundary = 'top'
    variable = 'disp_x'
    velocity = 'vel'
    component = 0
    associated_subdomain = 0
  []
  [displace_y_top]
    type = INSADDisplaceBoundaryBC
    boundary = 'top'
    variable = 'disp_y'
    velocity = 'vel'
    component = 1
    associated_subdomain = 0
  []
  [displace_x_top_dummy]
    type = INSADDummyDisplaceBoundaryIntegratedBC
    boundary = 'top'
    variable = 'disp_x'
    velocity = 'vel'
    component = 0
  []
  [displace_y_top_dummy]
    type = INSADDummyDisplaceBoundaryIntegratedBC
    boundary = 'top'
    variable = 'disp_y'
    velocity = 'vel'
    component = 1
  []
[]

[FunctorMaterials]
  [steel]
    type = AriaLaserWeld304LStainlessSteelFunctorMaterial
    temperature = T
    beta = 1e7
  []
  [disp_vec_value_and_dot]
    type = ADGenericVectorFunctorMaterial
    prop_names = 'disp_vec'
    prop_values = 'disp_x disp_y 0'
  []
  [vel]
    type = ADGenericVectorFunctorMaterial
    prop_names = 'vel'
    prop_values = 'vel_x vel_y 0'
  []
  [h]
    type = INSFVEnthalpyMaterial
    rho = 'rho'
    temperature = T
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_mat_solver_type -mat_mffd_err'
    petsc_options_value = 'lu       NONZERO               strumpack                  1e-6'
  []
[]

[Executioner]
  type = Transient
  end_time = ${endtime}
  dtmin = 1e-8
  dtmax = ${timestep}
  petsc_options = '-snes_converged_reason -ksp_converged_reason -options_left'
  solve_type = 'PJFNK'
  line_search = 'none'
  nl_max_its = 12
  l_max_its = 100
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 7
    dt = ${timestep}
    linear_iteration_ratio = 1e6
    growth_factor = 1.1
  []
[]

[Outputs]
  exodus = true
[]

[Debug]
  show_var_residual_norms = true
[]
