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
# pebble_bed_porosity       = 0.5 # (//).
pebble_bed_r              = 0.79788 # 2 m2 free flow area 1m2 of real flow area (m).
pebble_bed_h              = 10.0 # (m).
# pebble_bed_free_flow_area = ${fparse pi*pebble_bed_r*pebble_bed_r} # (m2).
# pebble_bed_free_volume    = ${fparse pebble_bed_free_flow_area * pebble_bed_h} # (m2).
pebbles_diameter          = 0.06 #(m).

# Properties -------------------------------------------------------------------
# total_power               = 0.0 # (W).
inlet_T_fluid             = 523.0 # (K).
outlet_pressure           = 7.0e+06 # (Pa).

# pebble_bed_power_density  = ${fparse total_power/pebble_bed_free_volume} # (W/m3)
inlet_superficial_vel     = 0.05 # 0.05 superficial velocity 0.1 real velocity (m/s)
rho_initial=1.29
rho_u_in=${fparse -rho_initial * inlet_superficial_vel}

user_limiter='min_mod'

[GlobalParams]
  pebble_diameter = ${pebbles_diameter}
  fp = fp
  two_term_boundary_expansion = true
  limiter = ${user_limiter}
  porosity = 1.0
  secondary_fraction = '0.1'
[]

[Mesh]
  type = MeshGeneratorMesh

  block_id = ' 1 '
  block_name = 'pebble_bed'

  uniform_refine = 2

  [cartesian_mesh]
    type = CartesianMeshGenerator
    dim = 2
    dx = '${pebble_bed_r}'
    ix = 4

    dy = '${pebble_bed_h}'
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
    initial_condition = 1e-15
  []
  [rho_v]
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
[]

[Modules]
  [FluidProperties]
    [fp_helium]
      type = IdealGasFluidProperties
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
    fraction = '0.1'
  []
  [porosity]
    type = GenericConstantMaterial
    prop_names = 'porosity'
    prop_values = '1.0'
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
  solve_type = NEWTON
  petsc_options = ' -snes_converged_reason '
  line_search = none

  # Problem time parameters.
  start_time = 0.0

  end_time = 1e+6
  dtmin    = 1e-4
  dtmax    = 5e+4


  # Iterations parameters.
  l_max_its = 50
  l_tol     = 1e-3

  nl_max_its = 25
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-7

  # Automatic scaling
  # automatic_scaling = true
  # compute_scaling_once = false

  # Steady state detection.
  steady_state_detection = true
  steady_state_tolerance = 1e-11
  # steady_state_tolerance = 1e-12
  # steady_state_start_time = 1e+3

  # Time step control.
  [TimeStepper]
    type = IterationAdaptiveDT
    dt                 = 2.5e-04
    cutback_factor     = 0.25
    growth_factor      = 2.00
    optimal_iterations = 12
  []

  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'
[]

[Outputs]
  [out]
    type = Exodus
    execute_on = 'initial timestep_end'
  []
[]