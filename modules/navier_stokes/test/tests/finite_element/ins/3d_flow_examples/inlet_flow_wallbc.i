[GlobalParams]
  order = FIRST
  family = LAGRANGE

  rho = rho
  u = u
  v = v
  w = w
  pressure = p
  temperature = T
  eos = eos

  conservative_form = true
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  xmin = 0
  xmax = 0.1
  ymin = 0
  ymax = 0.1
  zmin = 0
  zmax = 0.2
  nx = 4
  ny = 4
  nz = 8
  elem_type = HEX8
[]

[FluidProperties]
  [./eos]
    type = SimpleFluidProperties
    density0 = 1000              # kg/m^3
    thermal_expansion = 0       # K^{-1}
    cp =  4186.0
    viscosity = 1e-3             # Pa-s, Re=rho*u*L/mu = 100*1*0.1/0.1 = 100
  [../]
[]

[Variables]
  # velocity
  [./u]
    scaling = 1.e-1
    initial_condition = 0.0
  [../]
  [./v]
    scaling = 1.e-1
    initial_condition = 0.0
  [../]
  [./w]
    scaling = 1.e-3
    initial_condition = 1.0
  [../]

  # Temperature
  [./T]
    scaling = 1.e-6
    initial_condition = 783
  [../]

  # Pressure
  [./p]
    scaling = 1.e-1
    initial_condition = 1.2e5
  [../]
[]

[AuxVariables]
  [./rho]
    initial_condition = 1000
  [../]
  [./v_mag]
  [../]
[]

[Materials]
  [./pm]
    type = MDFluidMaterial
  [../]
[]

[Kernels]
  # mass
  [./mass_time]
    type = FluidPressureTimeDerivative
    variable = p
  [../]
  [./mass]
    type = MDFluidMassKernel
    variable = p
  [../]

  # x-momentum
  [./x_momentum_time]
    type = FluidVelocityTimeDerivative
    variable = u
  [../]
  [./x_momentum_space]
    type = MDFluidMomentumKernel
    variable = u
    component = 0
  [../]

  # y-momentum
  [./y_momentum_time]
    type = FluidVelocityTimeDerivative
    variable = v
  [../]
  [./y_momentum_space]
    type = MDFluidMomentumKernel
    variable = v
    component = 1
  [../]

  # z-momentum
  [./z_momentum_time]
    type = FluidVelocityTimeDerivative
    variable = w
  [../]
  [./z_momentum_space]
    type = MDFluidMomentumKernel
    variable = w
    component = 2
  [../]

  # temperature
  [./temperature_time]
    type = FluidTemperatureTimeDerivative
    variable = T
  [../]
  [./temperature_space]
    type = MDFluidEnergyKernel
    variable = T
  [../]
[]

[AuxKernels]
  [./rho_aux]
    type = NSDensityAux
    variable = rho
  [../]
  [./v_aux]
    type = VectorMagnitudeAux
    variable = v_mag
    x = u
    y = v
    z = w
  [../]
[]


[BCs]
  # Mass eqn BCs: mass inlet (back), pressure outlet (front),
  #               zero fluxes on other surfaces (not needed explicitly)
  [./inlet_mass]
    type = MDFluidMassBC
    variable = p
    boundary = 'back'
    v_fn = 1
  [../]
  [./pressure_outlet]
    type = DirichletBC
    variable = p
    boundary = 'front'
    value = 1e5
  [../]

  # x-momentum (u) BCs: zero on 'back right left bottom top'
  #                     pressure outlet on 'front'
  [./u_wall]
    type = DirichletBC
    variable = u
    boundary = 'back right left bottom top'
    value = 0
  [../]
  [./pressure_outlet_mom_x]
    type = MDFluidMomentumBC
    variable = u
    boundary = 'front'
    p_fn = 1e5
    component = 0
  [../]

  # y-momentum (v) BCs: zero on 'back right left bottom top'
  #                     pressure outlet on 'front'
  [./v_wall]
    type = DirichletBC
    variable = v
    boundary = 'back right left bottom top'
    value = 0.0
  [../]
  [./pressure_outlet_mom_y]
    type = MDFluidMomentumBC
    variable = v
    boundary = 'front'
    p_fn = 1e5
    component = 1
  [../]

  # z-momentum (w) BCs: zero on 'right left bottom top'
  #                     w-inlet (DirichletBC) on 'back'
  #                     pressure outlet on 'front'
  [./w_in]
    type = DirichletBC
    variable = w
    boundary = 'back'
    value = 1.0
  [../]
  [./z_mom_wall_bc]
    type = FluidWallMomentumBC
    variable = w
    boundary = 'right left bottom top'
    component = 2
  [../]
  [./pressure_outlet_mom_z]
    type = MDFluidMomentumBC
    variable = w
    boundary = 'front'
    p_fn = 1e5
    component = 2
  [../]

  # Energy equation BCs: T_inlet (DirichletBC) on 'back'
  #                      pressure outlet on 'front'
  #                      zero fluxes on other surfaces (not needed explicitly)
  #
  # note: it is probably more appropriate to use an integral BC for inlet condition,
  #       however, we keep what was originally setup
  [./T_bottom]
    type = DirichletBC
    variable = T
    boundary = 'back'
    value = 753
  [../]
  [./pressure_outlet_energy]
    type = MDFluidEnergyBC
    variable = T
    boundary = 'front'
  [../]
[]

[Preconditioning]
  [./SMP_PJFNK]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    dt = 1.e-1
  []
  petsc_options_iname = '-pc_type -ksp_gmres_restart'
  petsc_options_value = 'lu 100'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  end_time = 1e3
[]

[Outputs]
  perf_graph = true
  exodus = true
  print_linear_residuals = false
[]
