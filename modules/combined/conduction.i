[Mesh]
  [mesh]
    type = FileMeshGenerator
    file = 'mesh_conduction_in.e'
  []
[]

[Variables]
  [T]
    initial_condition = 1000
  []
[]

[AuxVariables]
  [power_density]
    family = LAGRANGE
    order = FIRST
  []
[]

[ICs]
  [power]
    type = ConstantIC
    variable = power_density
    value = 7675296.2525
    block = 'FUEL_11 FUEL_12 FUEL_13 FUEL_21 FUEL_22 FUEL_23
             FUEL_TRI_11 FUEL_TRI_12 FUEL_TRI_13 FUEL_TRI_21 FUEL_TRI_22 FUEL_TRI_23'
  []
[]

[Kernels]
  [conduction]
    type = HeatConduction
    variable = T
  []
  [heat]
    type = CoupledForce
    variable = T
    v = power_density
  []
[]

[Materials]
  [fuel]
    type = HeatConductionMaterial
    block = 'FUEL_11 FUEL_12 FUEL_13 FUEL_21 FUEL_22 FUEL_23
             FUEL_TRI_11 FUEL_TRI_12 FUEL_TRI_13 FUEL_TRI_21 FUEL_TRI_22 FUEL_TRI_23'
    thermal_conductivity = 1000
  []
  [monolith]
    type = HeatConductionMaterial
    block = 'MONOLITH MONOLITH_TRI'
    thermal_conductivity = 1830
  []
  [reflector]
    type = HeatConductionMaterial
    block = 'REFLECTOR REFLECTOR_TRI'
    thermal_conductivity = 200
  []
  [gap]
    type = HeatConductionMaterial
    block = 'GAP'
    thermal_conductivity = 0.08
  []
[]

[BCs]
  [htpipe]
    type = CoupledConvectiveHeatFluxBC
    variable = T
    boundary = 'htpipe'
    T_infinity = 800
    htc = 372
  []
  [convection]
    type = ConvectiveHeatFluxBC
    variable = T
    boundary = 'bottom top outer'
    T_infinity = 300
    heat_transfer_coefficient = 30
  []
  [radiation]
    type = RadiativeHeatFluxBC
    variable = T
    boundary = 'bottom top outer'
    Tinfinity = 300
    boundary_emissivity = 0.7
  []
[]

[Executioner]
  type = Steady
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 300'
  line_search = 'none'
  l_tol = 1e-02
  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-8
[]

[Outputs]
  perf_graph = true
  [console]
    type = Console
    output_file = true
  []
[]