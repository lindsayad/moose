# ==============================================================================
# Model description
# ------------------------------------------------------------------------------
# Idaho Falls, INL, June 15, 2021
# Author(s):Stefano Terlizzi
# ==============================================================================
# - iloop: Straight channel model.
# - mm: mass, momentum, energy.
# - rz: rz cohordinate system.
# ==============================================================================
# MODEL PARAMETERS
# ==============================================================================

# Geometry ---------------------------------------------------------------------
multiplier = 1e2
pebble_bed_porosity       = 0.5 # (//).
pebble_bed_r              = ${fparse 0.79788 * multiplier} # 2 m2 free flow area 1m2 of real flow area (m).
pebble_bed_h              = ${fparse 10.0 * multiplier} # (m).
# pebble_bed_free_flow_area = ${fparse pi*pebble_bed_r*pebble_bed_r} # (m2).
# pebble_bed_free_volume    = ${fparse pebble_bed_free_flow_area * pebble_bed_h} # (m2).
# pebbles_diameter          = 0.06 #(m).

# Properties -------------------------------------------------------------------
# total_power               = 0.0 # (W).
inlet_T_fluid             = 523.0 # (K).
outlet_pressure           = 7.0e+06 # (Pa).

# pebble_bed_power_density  = ${fparse total_power/pebble_bed_free_volume} # (W/m3)
inlet_superficial_vel     = ${fparse 0.05 * multiplier} # 0.05 superficial velocity 0.1 real velocity (m/s)
rho_initial = 46.6832

# T=273
# p_initial=1.01e5
# v_in=1
# gamma=1.4
# e_initial=${fparse p_initial / (gamma - 1) / rho_initial}
# et_initial=${e_initial}
# rho_et_initial=${fparse rho_initial * et_initial}
# rho_u_initial=${fparse rho_initial * v_in}

rho_u_in=${fparse -inlet_superficial_vel * rho_initial}

user_limiter='min_mod'
# user_limiter = 'upwind'

[GlobalParams]
  # pebble_diameter = ${pebbles_diameter}
  # acceleration = ' 0.00 -9.81 0.00 ' # Gravity acceleration (m/s2).
  fp = fp
  two_term_boundary_expansion = true
  limiter = ${user_limiter}
  porosity = ${pebble_bed_porosity}
[]

[Mesh]
  type = MeshGeneratorMesh

  block_id = ' 1 '
  block_name = 'pebble_bed'

  uniform_refine = 1

  [cartesian_mesh]
    type = CartesianMeshGenerator
    dim = 2
    dx = '${fparse pebble_bed_r}'
    ix = 4

    dy = '${fparse pebble_bed_h}'
    iy = '50'

    subdomain_id = ' 1 '
  []
[]

[Problem]
  coord_type = RZ
  material_coverage_check = false
  fv_bcs_integrity_check = false
[]

[Variables]
  [pressure]
    type = MooseVariableFVReal
    initial_condition = ${outlet_pressure}
  []
  [rho_u]
    type = MooseVariableFVReal
    initial_condition = 0.0
  []
  [rho_v]
    type = MooseVariableFVReal
    initial_condition = ${rho_u_in}
    # initial_condition = 1e-15
  []
  [Mass_Fraction]
    type = MooseVariableFVReal
    initial_condition = 1e-15
  []

[]

[FVKernels]
  # Mass conservation equation
  [mass_time]
    type = FVMatPropTimeKernel
    variable = pressure
    mat_prop_time_derivative = 'dsuperficial_rho_dt'
  []
  [mass_advection]
    type = GasMixPCNSFVKT
    variable = pressure
    eqn = "mass"
  []
  # Momentum conservation, x component
  [momentum_time_x]
    type = FVMatPropTimeKernel
    variable = rho_u
    mat_prop_time_derivative = 'dsuperficial_rhou_dt'
  []
  [momentum_advection_and_pressure_x]
    type = GasMixPCNSFVKT
    variable = rho_u
    eqn = "momentum"
    momentum_component = 'x'
  []
  [rz_pressure]
    type = PCNSFVMomentumPressureRZ
    variable = rho_u
  []
  [momentum_gravity_x]
    type = NSFVMomentumGravity
    variable = rho_u
    momentum_component = 'x'
    gravity = ' 0.00 -9.81 0.00 '
  []

  # Momentum conservation, y component

  [momentum_time_y]
    type = FVMatPropTimeKernel
    variable = rho_v
    mat_prop_time_derivative = 'dsuperficial_rhov_dt'
  []
  [momentum_advection_and_pressure_y]
    type = GasMixPCNSFVKT
    variable = rho_v
    eqn = "momentum"
    momentum_component = 'y'
  []
  [momentum_gravity_y]
    type = NSFVMomentumGravity
    variable = rho_v
    momentum_component = 'y'
    gravity = ' 0.00 -9.81 0.00 '
  []

  # Mass fraction advection
  [mass_frac_time]
    type = FVMatPropTimeKernel
    variable = Mass_Fraction
    mat_prop_time_derivative = 'drho_f_dt'
  []
  [MF_Advection]
    type = GasMixPCNSFVKT
    variable = Mass_Fraction
    eqn = "scalar"
  []
[]

[AuxVariables]
  [rho]
    type = MooseVariableFVReal
    initial_condition = ${rho_initial}
  []
  [T_fluid]
    type = MooseVariableFVReal
    initial_condition = ${inlet_T_fluid}
  []
  [vel_y]
    type = MooseVariableFVReal
  []
[]

[AuxKernels]
  [rho]
    type = ADMaterialRealAux
    property = rho
    variable = rho
  []
  [vel_y]
    type = ADMaterialRealAux
    property = vel_y
    variable = vel_y
    execute_on = 'timestep_end'
  []
[]

[FVBCs]
  #Reactor Inlet
  [pressure_inlet]
    type = GasMixPCNSFVStrongBC
    velocity_function_includes_rho = true
    superficial_velocity = 'rho_u_in_1'
    boundary = 'top'
    variable = pressure
    eqn = 'mass'
    T_fluid = ${inlet_T_fluid}
  []
  [superficial_vel_x_inlet]
    type = GasMixPCNSFVStrongBC
    velocity_function_includes_rho = true
    boundary = 'top'
    variable = rho_u
    superficial_velocity = 'rho_u_in_1'
    eqn = 'momentum'
    momentum_component = 'x'
    T_fluid = ${inlet_T_fluid}
  []
  [superficial_vel_y_inlet]
    type = GasMixPCNSFVStrongBC
    velocity_function_includes_rho = true
    boundary = 'top'
    variable = rho_v
    superficial_velocity = 'rho_u_in_1'
    eqn = 'momentum'
    momentum_component = 'y'
    T_fluid = ${inlet_T_fluid}
  []

  [Mass_Fraction_Inlet]
    type = GasMixPCNSFVStrongBC
    velocity_function_includes_rho = true
    boundary = 'top'
    variable = Mass_Fraction
    superficial_velocity = 'rho_u_in_1'
    eqn = 'scalar'
    scalar = 1
    T_fluid = ${inlet_T_fluid}
  []

  # Reactor Outlet Flow.
  [pressure_outlet]
    type = GasMixPCNSFVStrongBC
    boundary = 'bottom'
    variable = pressure
    p = ${outlet_pressure}
    eqn = 'mass'
  []

  [superficial_vel_x_outlet]
    type = GasMixPCNSFVStrongBC
    boundary = 'bottom'
    variable = rho_u
    p = ${outlet_pressure}
    eqn = 'momentum'
    momentum_component = 'x'
  []
  [superficial_vel_y_outlet]
    type = GasMixPCNSFVStrongBC
    boundary = 'bottom'
    variable = rho_v
    p = ${outlet_pressure}
    eqn = 'momentum'
    momentum_component = 'y'
  []
  [Mass_Fraction_Outler]
    type = GasMixPCNSFVStrongBC
    boundary = 'bottom'
    variable = Mass_Fraction
    p = ${outlet_pressure}
    eqn = 'scalar'
  []

  # Fluid solid walls.
  [superficial_vel_x_vertical_walls]
    type = PCNSFVImplicitMomentumPressureBC
    boundary = 'left right'
    variable = rho_u
    momentum_component = 'x'
  []
  [superficial_vel_y_vertical_walls_implicit]
    type = PCNSFVImplicitMomentumPressureBC
    boundary = 'left right'
    variable = rho_v
    momentum_component = 'y'
  []
  # Help gradient reconstruction
  [sup_mom_x_inlet_walls]
    type = FVDirichletBC
    boundary = 'top left right'
    variable = rho_u
    value = 0
  []
  [sup_mom_y_inlet]
    type = FVDirichletBC
    boundary = 'top'
    variable = rho_v
    value = ${rho_u_in}
  []
  [pressure_outlet_diri]
    type = FVDirichletBC
    variable = pressure
    boundary = 'bottom'
    value = ${outlet_pressure}
  []
  [MF_inlet]
    type = FVDirichletBC
    variable = Mass_Fraction
    boundary = 'top'
    value = 1
  []
[]

[Modules]
  [FluidProperties]
    [fp_helium]
      type = HeliumFluidProperties
    []
    [fp_air]
      type = IdealGasFluidProperties
    []
    [fp]
      type = GasMixPHFluidProperties
      fp_primary = fp_helium
      fp_secondary = 'fp_air'
    []
  []
[]

[Materials]
  [var_mat]
    type = GasMixPorousMixedVarMaterial
    p = pressure
    T_fluid = T_fluid
    superficial_rhou = rho_u
    superficial_rhov = rho_v
    secondary_fraction = Mass_Fraction
  []
  [porosity]
    type = GenericConstantMaterial
    prop_names = 'porosity'
    prop_values = '${pebble_bed_porosity}'
  []
[]

[Functions]
  [rho_u_in_1]
    type = ParsedVectorFunction
    value_x = '0'
    value_y = '${rho_u_in}'
  []
[]

[Executioner]
  type = Transient # Pseudo transient to reach steady state.

  # Problem time parameters.
  start_time = 0.0

  end_time = 200
  dtmax    = 1

  # Iterations parameters.
  l_max_its = 50
  l_tol     = 1e-5

  # Steady state detection.
  steady_state_detection = true
  steady_state_tolerance = 1e-11

  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'

  line_search = 'none'
  nl_max_its = 20
  nl_abs_tol = 1e-4
  solve_type = NEWTON
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 10
    dt = 1e-3
  []

  # num_steps = 1000
  # petsc_options = '-ksp_monitor'
  # dt = 1e-5
  # [TimeIntegrator]
  #   type = ExplicitSSPRungeKutta
  #   order = 2
    # type = ActuallyExplicitEuler
    # use_constant_mass = true
  # []
[]

[Outputs]
  checkpoint = true
  print_linear_residuals = true
  [out]
    type = Exodus
    execute_on = 'initial timestep_end'
  []
[]
