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
pebble_bed_porosity       = 0.5 # (//).
pebble_bed_r              = 0.79788 # 2 m2 free flow area 1m2 of real flow area (m).
pebble_bed_h              = 10.0 # (m).
pebble_bed_free_flow_area = ${fparse pi*pebble_bed_r*pebble_bed_r} # (m2).
pebble_bed_free_volume    = ${fparse pebble_bed_free_flow_area * pebble_bed_h} # (m2).
pebbles_diameter          = 0.06 #(m).

# Properties -------------------------------------------------------------------
total_power               = 0.0 # (W).
inlet_T_fluid             = 523.0 # (K).
outlet_pressure           = 7.0e+06 # (Pa).
pebble_bed_power_density  = ${fparse total_power/pebble_bed_free_volume} # (W/m3)
inlet_superficial_vel     = 0.05 # 0.05 superficial velocity 0.1 real velocity (m/s)
rho_initial=1.29

T=273
p_initial=1.01e6
v_in=1
gamma=1.4
e_initial=${fparse p_initial / (gamma - 1) / rho_initial}
et_initial=${e_initial}
rho_et_initial=${fparse rho_initial * et_initial}
rho_v_initial=${fparse rho_initial * v_in}

[GlobalParams]
  pebble_diameter = ${pebbles_diameter}
  acceleration = ' 0.00 -9.81 0.00 ' # Gravity acceleration (m/s2).
  fp = fp
  porosity = 1.0
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
  kernel_coverage_check = false
  material_coverage_check = false
[]

[Variables]
  [dummy]
    type = MooseVariableFVReal
    initial_condition = 0.0
  []
[]

[FVKernels]
  # Mass conservation equation
  [mass_time]
    type = FVPorosityTimeDerivative
    variable = dummy
  []
[]  

[ICs]
[]

[AuxVariables]
  [p]
    type = MooseVariableFVReal
  []
  [rho]
    type = MooseVariableFVReal
  []
  [rhou]
    type = MooseVariableFVReal
  []
  [rhov]
    type = MooseVariableFVReal
  []
  [T_fluid]
    type = MooseVariableFVReal
  []
  [p_mat]
    type = MooseVariableFVReal
  []
  [rho_mat]
    type = MooseVariableFVReal
  []
  [rhou_mat]
    type = MooseVariableFVReal
  []
  [rhov_mat]
    type = MooseVariableFVReal
  []
  [temperature_mat]
    type = MooseVariableFVReal
  []    
[]

[AuxKernels]
  [p]
    type = FunctionAux
    function = '${p_initial}'
    variable = p
    scaling = 1e-5
  []  
  [rho]
    type = FunctionAux
    function = '${rho_initial}'
    variable = rho
  []
  [rhou]
    type = FunctionAux
    function = 'x'
    variable = rhou
  []
  [rho_v]
    type = FunctionAux
    function = '${rho_v_initial}'
    variable = rhov
  []
  [T_fluid]
    type = FunctionAux
    function = '273.15'
    variable = T_fluid
  []
  [p_mat]
    type = ADMaterialRealAux
    property = p
    variable = p_mat
  []    
  [rho_mat]
    type = ADMaterialRealAux
    property=rho
    variable = rho_mat
  []
  [rhou_mat]
    type = ADMaterialRealAux
    property = rhou
    variable = rhou_mat
  []
  [rho_v_mat]
    type = ADMaterialRealAux
    property = rhov
    variable = rhov_mat
  []
  [temperature_mat]
    type = ADMaterialRealAux
    property = T_fluid
    variable = temperature_mat
  []     
[]

[FVBCs]
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
    p = p
    T_fluid = T_fluid
    superficial_rhou = rhou
    superficial_rhov = rhov
  []
  [porosity]
    type = GenericConstantMaterial
    prop_names = 'porosity'
    prop_values = '1.0'
  []
[]

[Postprocessors]
  [rho_difference]
    type = ElementL2Difference
    variable = rho
    other_variable = rho_mat
  []
  [rhou_difference]
    type = ElementL2Difference
    variable = rhou
    other_variable = rhou_mat
  []
  [rhov_difference]
    type = ElementL2Difference
    variable = rhov
    other_variable = rhov_mat
  []
  [p_difference]
    type = ElementL2Difference
    variable = p
    other_variable = p_mat
  []
  [temperature_difference]
    type = ElementL2Difference
    variable = T_fluid
    other_variable = temperature_mat
  []
[]

[Executioner]
 # type = Transient # Pseudo transient to reach steady state.
 # solve_type = NEWTON
 # petsc_options = ' -snes_converged_reason '
 # line_search = none

 # Problem time parameters.
  start_time = 0.0
  end_time = 10
  dt= 1
  type = Transient # Pseudo transient to reach steady state.
  solve_type = NEWTON
  petsc_options = ' -snes_converged_reason '
  line_search = none
[]  

[Outputs]
  [out]
    type = Exodus
    execute_on = 'initial timestep_end'
  []
[]

