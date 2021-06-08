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
# pebbles_diameter          = 0.06 #(m).

# Properties -------------------------------------------------------------------
# total_power               = 0.0 # (W).
inlet_T_fluid             = 523.0 # (K).
outlet_pressure           = 7.0e+06 # (Pa).

# pebble_bed_power_density  = ${fparse total_power/pebble_bed_free_volume} # (W/m3)
# inlet_superficial_vel     = 0.05 # 0.05 superficial velocity 0.1 real velocity (m/s)
inlet_superficial_vel = 500
rho_initial=1.29
rho_u_in=${fparse -rho_initial * inlet_superficial_vel}

user_limiter='upwind'

[GlobalParams]
  # pebble_diameter = ${pebbles_diameter}
  # fpressure = fp
  two_term_boundary_expansion = true
  limiter = ${user_limiter}
  porosity = 1.0
  # secondary_fraction = '0.5'
  fp = fp
[]

[Mesh]
  type = MeshGeneratorMesh

  block_id = ' 1 '
  block_name = 'pebble_bed'

  uniform_refine = 0

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
  [mass_advection]
    type = PCNSFVKT
    variable = pressure
    eqn = "mass"
  []
  # Momentum conservation, x component
  [momentum_advection_and_pressure_x]
    type = PCNSFVKT
    variable = rho_u
    eqn = "momentum"
    momentum_component = 'x'
  []
  [rz_pressure]
    type = PCNSFVMomentumPressureRZ
    variable = rho_u
  []

  # Momentum conservation, y component

  [momentum_advection_and_pressure_y]
    type = PCNSFVKT
    variable = rho_v
    eqn = "momentum"
    momentum_component = 'y'
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
    type = PCNSFVStrongBC
    velocity_function_includes_rho = true
    superficial_velocity = 'rho_u_in_1'
    boundary = 'top'
    variable = pressure
    eqn = 'mass'
    T_fluid = ${inlet_T_fluid}
  []
  [superficial_vel_x_inlet]
    type = PCNSFVStrongBC
    velocity_function_includes_rho = true
    boundary = 'top'
    variable = rho_u
    superficial_velocity = 'rho_u_in_1'
    eqn = 'momentum'
    momentum_component = 'x'
    T_fluid = ${inlet_T_fluid}
  []
  [superficial_vel_y_inlet]
    type = PCNSFVStrongBC
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
    type = PCNSFVStrongBC
    boundary = 'bottom'
    variable = pressure
    pressure = ${outlet_pressure}
    eqn = 'mass'
  []

  [superficial_vel_x_outlet]
    type = PCNSFVStrongBC
    boundary = 'bottom'
    variable = rho_u
    pressure = ${outlet_pressure}
    eqn = 'momentum'
    momentum_component = 'x'
  []
  [superficial_vel_y_outlet]
    type = PCNSFVStrongBC
    boundary = 'bottom'
    variable = rho_v
    pressure = ${outlet_pressure}
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
    [fp]
      type = IdealGasFluidProperties
    []
  []
[]

[Materials]
  [var_mat]
    type = PorousMixedVarMaterial
    p = pressure
    T_fluid = T_fluid
    superficial_rhou = rho_u
    superficial_rhov = rho_v
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

[Postprocessors]
  [average_density]
    type = ADElementAverageMaterialProperty
    mat_prop = rho
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options = ' -snes_converged_reason '

  nl_rel_tol = 1e-8

  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'
[]

[Outputs]
  exodus = true
[]
