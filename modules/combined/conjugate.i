# Units: specific_heat_capacity--cp--J/(kg.K); density--rho--kg/(cm^3);
# dynamic_viscosity--mu--kg/(cm.s); thermal_conductivity--k--W/(cm.K);
# pressure--kg/(cm.s^2); force--kg.cm/s^2

outlet_pressure = 0
inlet_velocity = 150 # cm/s
ini_temp = 300 # K
heat_transfer_coefficient = 9 # W/(cm2.K)
g = -981 # cm/s2
alpha_fluid = 2e-4 # thermal expansion coefficient of fluid used in INSADBoussinesqBodyForce

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  file = 'solid_fluid.msh'
[]

[Variables]
  [./T]
    initial_condition = ${ini_temp}
  [../]
  [./velocity]
    family = LAGRANGE_VEC
    block = 'fluid'
  [../]
  [./p]
    order = FIRST
    family = LAGRANGE
    block = 'fluid'
  [../]
[]

[AuxVariables]
  [./heat]
    family = LAGRANGE
    order = FIRST
    block = 'solid'
  [../]
[]

[ICs]
  [./initial_velocity]
    type = VectorConstantIC
    x_value = 0
    y_value = 0
    z_value = ${inlet_velocity}
    variable = velocity
  [../]
  [./initial_p]
    type = FunctionIC
    variable = p
    function = ini_p
  [../]
[]

[Kernels]
  [./fluid_mass]
    type = INSADMass
    variable = p
  [../]
  [./fluid_mass_pspg]
    type = INSADMassPSPG
    variable = p
  [../]

  [./fluid_momentum_time]
    type = INSADMomentumTimeDerivative
    variable = velocity
  [../]
  [./fluid_momentum_convection]
    type = INSADMomentumAdvection
    variable = velocity
  [../]
  [./fluid_momentum_viscous]
    type = INSADMomentumViscous
    variable = velocity
  [../]
  [./fluid_momentum_pressure]
    type = INSADMomentumPressure
    variable = velocity
    pressure = p
    integrate_p_by_parts = true
  [../]
  [./fluid_momentum_gravity]
    type = INSADGravityForce
    variable = velocity
    gravity = '0 0 ${g}'
  [../]
  [./fluid_momentum_buoyancy]
    type = INSADBoussinesqBodyForce
    variable = velocity
    gravity = '0 0 ${g}'
    alpha_name = 'alpha_fluid'
    ref_temp = 'T_ref'
	temperature = T
  [../]
  [./fluid_momentum_supg]
    type = INSADMomentumSUPG
    variable = velocity
    velocity = velocity
  [../]

  [./fluid_temperature_time]
    type = INSADHeatConductionTimeDerivative
    variable = T
    block = 'fluid'
  [../]
  [./fluid_temperature_conduction]
    type = ADHeatConduction
    variable = T
    block = 'fluid'
    thermal_conductivity = 'k'
  [../]
  [./fluid_temperature_advection]
    type = INSADEnergyAdvection
    variable = T
    block = 'fluid'
  [../]
  [./fluid_temperature_supg]
    type = INSADEnergySUPG
    variable = T
    velocity = velocity
    block = 'fluid'
  [../]

  [./solid_temperature_time]
    type = ADHeatConductionTimeDerivative
    variable = T
    block = 'solid'
    density_name = 'rho'
    specific_heat = 'cp'
  [../]
  [./solid_temperature_conduction]
    type = ADHeatConduction
    variable = T
    block = 'solid'
    thermal_conductivity = 'k'
  [../]
  [./heat_source]
    type = ADCoupledForce
    variable = T
    block = 'solid'
    v = heat
  [../]
[]

[AuxKernels]
  [./heat_aux]
    type = FunctionAux
    variable = heat
    function = heat_function
    block = 'solid'
  [../]
[]

[Modules/TensorMechanics/Master]
  [all]
    add_variables = true
    strain = FINITE
    automatic_eigenstrain_names = true
    generate_output = 'vonmises_stress stress_xx stress_yy stress_zz strain_xx strain_yy strain_zz'
    block = 'solid'
    material_output_order = FIRST
  []
[]

[InterfaceKernels]
  [./convection_heat_transfer]
    type = ConjugateHeatTransfer
    variable = T
    T_fluid = T
    neighbor_var = 'T'
    boundary = 'solid_wall'
    htc = 'htc'
  [../]
[]

[BCs]
  [./pin1_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'pin1'
    value = 0
  [../]
  [./pin1_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'pin1'
    value = 0
  [../]
  [./pin2_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'pin2'
    value = 0
  [../]
  [./solid_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'solid_bottom solid_top'
    value = 0
  [../]

  [./inlet_T_fluid]
    type = DirichletBC
    variable = T
    boundary = 'fluid_bottom'
    value = ${ini_temp}
  [../]
  [./inlet_velocity_fluid]
    type = VectorFunctionDirichletBC
    variable = velocity
    boundary = 'fluid_bottom'
    function_z = ${inlet_velocity}
  [../]
  [./outlet_p_fluid]
    type = DirichletBC
    variable = p
    boundary = 'fluid_top'
    value = ${outlet_pressure}
  [../]

  [./no_slip]
    type = VectorFunctionDirichletBC
    variable = velocity
    boundary = 'solid_wall'
  [../]

  [./Periodic]
    [./x]
      variable = velocity
      primary = 'fluid_wall1'
      secondary = 'fluid_wall2'
      translation = '1.6 0 0'
    [../]
    [./y]
      variable = velocity
      primary = 'fluid_wall3'
      secondary = 'fluid_wall4'
      translation = '0 1.6 0'
    [../]
  [../]
[]

[Materials]
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2e7 # N cm-2
    poissons_ratio = 0.32
    block = 'solid'
  []
  [thermal_expansion]
    type = ComputeThermalExpansionEigenstrain
    temperature = T
    thermal_expansion_coeff = 1.0358e-5
    stress_free_temperature = 273
    eigenstrain_name = thermal_expansion
    block = 'solid'
  []
  [stress]
    type = ComputeFiniteStrainElasticStress
    block = 'solid'
  []

  [./rho_solid]
    type = ADParsedMaterial
    f_name = rho
    function = '0.0110876 * pow(9.9672e-1 + 1.179e-5 * T - 2.429e-9 * pow(T,2) + 1.219e-12 * pow(T,3),-3)'
    args = 'T'
    block = 'solid'
  [../]
  [./cp_solid]
    type = ADParsedMaterial
    f_name = cp
    function = '0.76 * ((302.27 * pow((548.68/T),2) * exp(548.68 / T)) / pow((exp(548.68 / T) - 1),2) + 2 * 8.463e-3 * T + 8.741e7 * 18531.7 * exp(-18531.7 / T) / pow(T,2)) + 0.24 * ((322.49 * pow((587.41/T),2) * exp(587.41 / T)) / pow((exp(587.41 / T) - 1),2) + 2 * 1.4679e-2 * T)'
    args = 'T'
    block = 'solid'
  [../]
  [./k_solid]
    type = ADParsedMaterial
    f_name = k
    function = '1.158/(7.5408 + 17.692 * (T / 1000) + 3.6142 * pow((T/1000),2)) + 74.105 * pow((T / 1000),-2.5) * exp(-16.35 / (T / 1000))'
    args = 'T'
    block = 'solid'
  [../]

  [./rho_fluid]
    type = ADParsedMaterial
    f_name = rho
    function = '(11096 - 1.3236 * T) * 1e-6'
    args = 'T'
    block = 'fluid'
  [../]
  [./cp_fluid]
    type = ADParsedMaterial
    f_name = cp
    function = '159 - 2.72e-2 * T + 7.12e-6 * pow(T,2)'
    args = 'T'
    block = 'fluid'
  [../]
  [./k_fluid]
    type = ADParsedMaterial
    f_name = k
    function = '(3.61 + 1.517e-2 * T - 1.741e-6 * pow(T,2)) * 1e-2'
    args = 'T'
    block = 'fluid'
  [../]
  [./mu_fluid]
    type = ADParsedMaterial
    f_name = mu
    function = '4.94e-6 * exp(754.1/T)'
    args = 'T'
    block = 'fluid'
  [../]
  [./buoyancy_thermal_expansion_coefficient_fluid]
    type = ADGenericConstantMaterial
    prop_names = 'alpha_fluid'
    prop_values = '${alpha_fluid}'
	block = 'fluid'
  [../]
  [./buoyancy_reference_temperature_fluid]
    type = GenericConstantMaterial
    prop_names = 'T_ref'
    prop_values = '${ini_temp}'
	block = 'fluid'
  [../]

  [./ins_mat_fluid]
    type = INSADStabilized3Eqn
    velocity = velocity
    pressure = p
    temperature = T
    block = 'fluid'
  [../]

  [./htc]
    type = ADGenericFunctionMaterial
    prop_names = htc
    prop_values = htc_function
  [../]
[]

[Functions]
  [./htc_function]
    type = ParsedFunction
    value = ${heat_transfer_coefficient}
  [../]
  [./ini_p]
    type = ParsedFunction
    value = '0.01*981*(10-z)'
  [../]
  [./heat_function]
    type = ParsedFunction
    value = (100+200*sin((z/5)*(pi/2)))
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    solve_type = 'NEWTON'
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 100
  dt = 0.01
  #end_time = 100

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10

  solve_type = 'NEWTON'
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'
  line_search = 'none'
   # petsc_options_iname = '-snes_type'
  # petsc_options_value = 'test'

  nl_max_its = 30
  l_max_its = 100
  automatic_scaling = true
[]

[Postprocessors]
  [./max_T_solid]
    type = ElementExtremeValue
    variable = T
    value_type = max
    block = 'solid'
  [../]
  [./max_T_fluid]
    type = ElementExtremeValue
    variable = T
    value_type = max
    block = 'fluid'
  [../]
  [./min_T_solid]
    type = ElementExtremeValue
    variable = T
    value_type = min
    block = 'solid'
  [../]
  [./min_T_fluid]
    type = ElementExtremeValue
    variable = T
    value_type = min
    block = 'fluid'
  [../]
  [./average_T_solid]
    type = ElementAverageValue
    variable = T
    block = 'solid'
  [../]
  [./average_T_fluid]
    type = ElementAverageValue
    variable = T
    block = 'fluid'
  [../]
[]

[Outputs]
  perf_graph = true
  print_linear_residuals = true
  [./exodus]
    type = Exodus
    file_base = 'couple'
    execute_on = 'TIMESTEP_END'
  [../]
  [./csv]
    type = CSV
    file_base = 'couple'
    execute_on = 'TIMESTEP_END'
  [../]
[]
