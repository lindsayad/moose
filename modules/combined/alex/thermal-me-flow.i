# Units: specific_heat_capacity--cp--J/(kg.K); density--rho--kg/(cm^3);
# dynamic_viscosity--mu--kg/(cm.s); thermal_conductivity--k--W/(cm.K);
# pressure--kg/(cm.s^2); force--kg.cm/s^2

outlet_pressure = 0
inlet_velocity = 150 # cm/s
ini_temp = 593 # K
heat_transfer_coefficient = 9 # W/(cm2.K)
g = -981 # cm/s2
alpha_outer = 2e-4 # thermal expansion coefficient of outer used in INSADBoussinesqBodyForce

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  file = '3layers_3d_4parts_fine.msh'
[]

[Variables]
  [velocity]
    family = LAGRANGE_VEC
    order = FIRST
    block = 'outer'
  []
  [p]
    family = LAGRANGE
    order = FIRST
    block = 'outer'
  []
  [Tf]
    family = LAGRANGE
    order = FIRST
    block = 'outer'
  []
  [Ts]
    family = LAGRANGE
    order = FIRST
    block = 'inner mid'
  []
  [disp_x]
    family = LAGRANGE
    order = FIRST
    block = 'inner mid outer'
  []
  [disp_y]
    family = LAGRANGE
    order = FIRST
    block = 'inner mid outer'
  []
  [disp_z]
    family = LAGRANGE
    order = FIRST
    block = 'inner mid outer'
  []
[]

[AuxVariables]
  [power]
    family = MONOMIAL
    order = FIRST
    block = 'inner'
  []
[]

[ICs]
  [initial_velocity]
    type = VectorConstantIC
    variable = velocity
    x_value = 0
    y_value = 0
    z_value = ${inlet_velocity}
  []
  [initial_p]
    type = FunctionIC
    variable = p
    function = ini_p
  []
  [initial_Tf]
    type = ConstantIC
    variable = Tf
    value = ${ini_temp}
  []
  [initial_Ts]
    type = ConstantIC
    variable = Ts
    value = ${ini_temp}
  []
[]

[Kernels]
  [outer_mass]
    type = INSADMass
    variable = p
    use_displaced_mesh = true
  []
  [outer_mass_pspg]
    type = INSADMassPSPG
    variable = p
    use_displaced_mesh = true
  []
  [outer_momentum_time]
    type = INSADMomentumTimeDerivative
    variable = velocity
    use_displaced_mesh = true
  []
  [outer_momentum_convection]
    type = INSADMomentumAdvection
    variable = velocity
    use_displaced_mesh = true
  []
  [outer_momentum_viscous]
    type = INSADMomentumViscous
    variable = velocity
    use_displaced_mesh = true
  []
  [outer_momentum_pressure]
    type = INSADMomentumPressure
    variable = velocity
    pressure = p
    integrate_p_by_parts = true
    use_displaced_mesh = true
  []
  [outer_momentum_gravity]
    type = INSADGravityForce
    variable = velocity
    gravity = '0 0 ${g}'
    use_displaced_mesh = true
  []
  [outer_momentum_buoyancy]
    type = INSADBoussinesqBodyForce
    variable = velocity
    gravity = '0 0 ${g}'
    alpha_name = 'alpha_outer'
    ref_temp = 'T_ref'
    temperature = Tf
    use_displaced_mesh = true
  []
  [outer_momentum_supg]
    type = INSADMomentumSUPG
    variable = velocity
    velocity = velocity
    use_displaced_mesh = true
  []
  [outer_temperature_time]
    type = INSADHeatConductionTimeDerivative
    variable = Tf
    use_displaced_mesh = true
  []
  [outer_temperature_conduction]
    type = ADHeatConduction
    variable = Tf
    thermal_conductivity = 'k'
    use_displaced_mesh = true
  []
  [outer_temperature_advection]
    type = INSADEnergyAdvection
    variable = Tf
    use_displaced_mesh = true
  []
  [outer_temperature_supg]
    type = INSADEnergySUPG
    variable = Tf
    velocity = velocity
    use_displaced_mesh = true
  []

  [mid_inner_temperature_time]
    type = ADHeatConductionTimeDerivative
    variable = Ts
    density_name = 'rho'
    specific_heat = 'cp'
    block = 'mid inner'
    use_displaced_mesh = true
  []
  [mid_inner_temperature_conduction]
    type = ADHeatConduction
    variable = Ts
    thermal_conductivity = 'k'
    block = 'mid inner'
    use_displaced_mesh = true
  []
  [heat_source]
    type = ADCoupledForce
    variable = Ts
    v = power
    block = 'inner'
    use_displaced_mesh = true
  []

  [disp_x_diff]
    type = Diffusion
    variable = disp_x
    block = outer
    use_displaced_mesh = true
  []
  [disp_y_diff]
    type = Diffusion
    variable = disp_y
    block = outer
    use_displaced_mesh = true
  []
  [disp_z_diff]
    type = Diffusion
    variable = disp_z
    block = outer
    use_displaced_mesh = true
  []
[]

[Modules/TensorMechanics/Master]
  displacements = 'disp_x disp_y disp_z'
  strain = FINITE
  material_output_order = FIRST
  generate_output = 'vonmises_stress stress_xx stress_yy stress_zz strain_xx strain_yy strain_zz'
  [mechanics]
    block = 'mid inner'
    temperature = Ts
    displacements = 'disp_x disp_y disp_z'
    automatic_eigenstrain_names = true
  []
[]

[InterfaceKernels]
  [convection_heat_transfer]
    type = ConjugateHeatTransfer
    variable = Tf
    T_fluid = Tf
    neighbor_var = 'Ts'
    boundary = 'mid_wall'
    htc = 'htc'
    use_displaced_mesh = true
  []
[]

[AuxKernels]
  [power_distribution_auxk]
    type = FunctionAux
    variable = power
    function = power_distribution_function
    block = 'inner'
    use_displaced_mesh = true
  []
[]

[BCs]
  [no_slip]
    type = VectorFunctionDirichletBC
    variable = velocity
    boundary = 'mid_wall'
    use_displaced_mesh = true
  []
  [inlet_velocity]
    type = VectorFunctionDirichletBC
    variable = velocity
    boundary = 'outer_bottom'
    function_z = ${inlet_velocity}
    use_displaced_mesh = true
  []
  [symmetry]
    type = ADVectorFunctionDirichletBC
    variable = velocity
    boundary = 'outer_wall'
    function_x = 0
    function_y = 0
    set_x_comp = true
    set_y_comp = true
    set_z_comp = false
    use_displaced_mesh = true
  []
  [outlet_p]
    type = DirichletBC
    variable = p
    boundary = 'outer_top'
    value = ${outlet_pressure}
    use_displaced_mesh = true
  []
  [inlet_T]
    type = DirichletBC
    variable = Tf
    boundary = 'outer_bottom'
    value = ${ini_temp}
    use_displaced_mesh = true
  []

  [pin1_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'pin1'
    value = 0
    use_displaced_mesh = true
  []
  [pin1_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'pin1'
    value = 0
    use_displaced_mesh = true
  []
  [pin2_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'pin2'
    value = 0
    use_displaced_mesh = true
  []
  [pin2_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'pin2'
    value = 0
    use_displaced_mesh = true
  []
  [pin3_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'pin3'
    value = 0
    use_displaced_mesh = true
  []
  [axial_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'inner_bottom mid_bottom outer_bottom'
    value = 0
    use_displaced_mesh = true
  []
  [radial_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'outer_wall outer_bottom'
    value = 0
    use_displaced_mesh = true
  []
  [radial_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'outer_wall outer_bottom'
    value = 0
    use_displaced_mesh = true
  []
[]

[Materials]
  [rho_inner]
    type = ADParsedMaterial
    f_name = rho
    function = '0.0110876 * pow(9.9672e-1 + 1.179e-5 * Ts - 2.429e-9 * pow(Ts,2) + 1.219e-12 * pow(Ts,3),-3)'
    args = 'Ts'
    block = 'inner'
  []
  [cp_inner]
    type = ADParsedMaterial
    f_name = cp
    function = '0.76 * ((302.27 * pow((548.68 / Ts),2) * exp(548.68 / Ts)) / pow((exp(548.68 / Ts) - 1),2) + 2 * 8.463e-3 * Ts + 8.741e7 * 18531.7 * exp(-18531.7 / Ts) / pow(Ts,2)) + 0.24 * ((322.49 * pow((587.41/Ts),2) * exp(587.41 / Ts)) / pow((exp(587.41 / Ts) - 1),2) + 2 * 1.4679e-2 * Ts)'
    args = 'Ts'
    block = 'inner'
  []
  [k_inner]
    type = ADParsedMaterial
    f_name = k
    function = '1.158/(7.5408 + 17.692 * (Ts / 1000) + 3.6142 * pow((Ts/1000),2)) + 74.105 * pow((Ts / 1000),-2.5) * exp(-16.35 / (Ts / 1000))'
    args = 'Ts'
    block = 'inner'
  []

  [rho_mid]
    type = ADParsedMaterial
    f_name = rho
    function = '1e-6 * (7830.853 - 0.212046 * Ts - 1.011373e-4 * pow(Ts,2))'
    args = 'Ts'
    block = 'mid'
  []
  [cp_mid]
    type = ADParsedMaterial
    f_name = cp
    function = '5863.9 - 32.563 * Ts + 0.072564 * pow(Ts,2) - 7.045375e-5 * pow(Ts,3) + 2.585336e-8 * pow(Ts,4)'
    args = 'Ts'
    block = 'mid'
  []
  [k_mid]
    type = ADParsedMaterial
    f_name = k
    function = '0.3'
    args = 'Ts'
    block = 'mid'
  []

  [rho_outer]
    type = ADParsedMaterial
    f_name = rho
    function = '(11096 - 1.3236 * Tf) * 1e-6'
    args = 'Tf'
    block = 'outer'
  []
  [cp_outer]
    type = ADParsedMaterial
    f_name = cp
    function = '159 - 2.72e-2 * Tf + 7.12e-6 * pow(Tf,2)'
    args = 'Tf'
    block = 'outer'
  []
  [k_outer]
    type = ADParsedMaterial
    f_name = k
    function = '(3.61 + 1.517e-2 * Tf - 1.741e-6 * pow(Tf,2)) * 1e-2'
    args = 'Tf'
    block = 'outer'
  []
  [mu_outer]
    type = ADParsedMaterial
    f_name = mu
    function = '4.94e-6 * exp(754.1/Tf)'
    args = 'Tf'
    block = 'outer'
  []
  [buoyancy_thermal_expansion_coefficient_outer]
    type = ADGenericConstantMaterial
    prop_names = 'alpha_outer'
    prop_values = '${alpha_outer}'
    block = 'outer'
  []
  [buoyancy_reference_temperature_outer]
    type = GenericConstantMaterial
    prop_names = 'T_ref'
    prop_values = '${ini_temp}'
    block = 'outer'
  []

  [ins_mat_outer]
    type = INSADStabilized3Eqn
    velocity = velocity
    pressure = p
    temperature = Tf
    block = 'outer'
  []

  [htc]
    type = ADGenericFunctionMaterial
    prop_names = htc
    prop_values = htc_function
    use_displaced_mesh = true
  []

  [elasticity_inner_mid]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2e7
    poissons_ratio = 0.32
    block = 'inner mid'
    use_displaced_mesh = true
  []
  [thermal_expansion_inner_mid]
    type = ComputeThermalExpansionEigenstrain
    temperature = Ts
    thermal_expansion_coeff = 2e-5
    stress_free_temperature = 593
    eigenstrain_name = thermal_expansion
    block = 'inner mid'
    use_displaced_mesh = true
  []
  [stress_inner_mid]
    type = ComputeFiniteStrainElasticStress
    block = 'inner mid'
  []
[]

[Functions]
  [htc_function]
    type = ParsedFunction
    value = ${heat_transfer_coefficient}
  []
  [ini_p]
    type = ParsedFunction
    value = '0.010302 * 981 * (90 - z)'
  []
  [power_distribution_function]
    type = ParsedFunction
    value = '300 * sin(pi * z / 90)'
  []
[]

[Preconditioning]
  [fsp]
    type = FSP
    topsplit = 'up'
    # full = false
    [up]
      splitting = 'u p'
      splitting_type = schur
      petsc_options_iname = '-pc_fieldsplit_schur_fact_type  -pc_fieldsplit_schur_precondition'
      petsc_options_value = 'full selfp'
    []
    [u]
      vars = 'Ts Tf velocity p'
      petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_mat_solver_type -ksp_type'
      petsc_options_value = 'lu       NONZERO               superlu_dist               preonly'
    []
    [p]
      vars = 'disp_x disp_y disp_z'
      petsc_options_iname = '-pc_type -pc_hypre_type -ksp_type'
      petsc_options_value = 'hypre    boomeramg      preonly'
    []
  []
[]

[Executioner]
  type = Transient
  dt = 1
  num_steps = 3
  # end_time = 1
  steady_state_detection = true

  solve_type = 'NEWTON'
  # solve_type = 'PJFNK'
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor'

  # petsc_options_iname = '-pc_type -pc_factor_shift_type'
  # petsc_options_value = 'lu       NONZERO'

  # petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_levels -ksp_gmres_restart'
  # petsc_options_value = 'asm      ilu          3                     200               '

  # petsc_options_iname = '-pc_type -sub_pc_type'
  # petsc_options_value = 'asm      lu          '

  # petsc_options_iname = '-pc_type -sub_pc_type -pc_hypre_type'
  # petsc_options_value = 'asm      hypre        boomeramg     '

  # petsc_options_iname = '-pc_type -pc_hypre_type'
  # petsc_options_value = 'hypre euclid'

  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  # petsc_options_value = 'lu       mumps'

  # petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_interp_type -pc_hypre_boomeramg_truncfactor'
  # petsc_options_value = 'hypre    boomeramg      0.7                                    ext+i                           0.2                            '

  # petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_max_levels -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type -pc_hypre_boomeramg_P_max -pc_hypre_boomeramg_truncfactor'
  # petsc_options_value = 'hypre    boomeramg      0.7                                  3                          4                                  25                            HMIS                             ext+i                           2                         0.3                            '

  line_search = 'none'

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8

  nl_max_its = 30
  l_max_its = 100
  automatic_scaling = true
  compute_scaling_once = false
  off_diagonals_in_auto_scaling = true
[]

[Outputs]
  perf_graph = true
  print_linear_residuals = true
  [exodus]
    type = Exodus
    file_base = 'thermal-me'
    execute_on = 'TIMESTEP_END'
  []
  [csv]
    type = CSV
    file_base = 'thermal-me'
    execute_on = 'TIMESTEP_END'
  []
[]

[Postprocessors]
  [average_inner_Ts]
    type = ElementAverageValue
    variable = Ts
    block = 'inner'
    use_displaced_mesh = true
  []
  [average_mid_Ts]
    type = ElementAverageValue
    variable = Ts
    block = 'mid'
    use_displaced_mesh = true
  []
  [average_outer_Tf]
    type = ElementAverageValue
    variable = Tf
    block = 'outer'
    use_displaced_mesh = true
  []
  [max_inner_Ts]
    type = ElementExtremeValue
    variable = Ts
    value_type = max
    block = 'inner'
    use_displaced_mesh = true
  []
  [max_mid_Ts]
    type = ElementExtremeValue
    variable = Ts
    value_type = max
    block = 'mid'
    use_displaced_mesh = true
  []
  [max_outer_Tf]
    type = ElementExtremeValue
    variable = Tf
    value_type = max
    block = 'outer'
    use_displaced_mesh = true
  []
  [min_inner_Ts]
    type = ElementExtremeValue
    variable = Ts
    value_type = min
    block = 'inner'
    use_displaced_mesh = true
  []
  [min_mid_Ts]
    type = ElementExtremeValue
    variable = Ts
    value_type = min
    block = 'mid'
    use_displaced_mesh = true
  []
  [min_outer_Tf]
    type = ElementExtremeValue
    variable = Tf
    value_type = min
    block = 'outer'
    use_displaced_mesh = true
  []
[]

[Debug]
  show_var_residual_norms = true
  # show_actions=true
[]
