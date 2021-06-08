rho_initial=1.29
p_initial=1.01e5
T=273.15
gamma=1.4
e_initial=${fparse p_initial / (gamma - 1) / rho_initial}
et_initial=${e_initial}
rho_et_initial=${fparse rho_initial * et_initial}
v_in=.5
rho_v_initial=${fparse rho_initial * v_in}



[GlobalParams]
  fp = fp
[]

[Mesh]
  [cartesian]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 1
    nx = 2
    ymin = 0
    ymax = 10
    ny = 50
  []
[]

[Modules]
  [FluidProperties]
    [fp]
      type = IdealGasFluidProperties
    []
  []
[]

[Variables]
  [mass_frac]
    type = MooseVariableFVReal
  []
[]

[AuxVariables]
  [rho]
    type = MooseVariableFVReal
  []
  [U_x]
    type = MooseVariableFVReal
  []
  [U_y]
    type = MooseVariableFVReal
  []
  [pressure]
    type = MooseVariableFVReal
  []
  [temperature]
    type = MooseVariableFVReal
  []
  [courant]
    type = MooseVariableFVReal
  []
  [rho_u]
    type = MooseVariableFVReal
    initial_condition = 1e-15
  []
  [rho_v]
    type = MooseVariableFVReal
    initial_condition = ${rho_v_initial}
  []
  [rho_et]
    type = MooseVariableFVReal
    initial_condition = ${rho_et_initial}
  []
[]


[ICs]
  [mf_initial]
    type = FunctionIC
    variable = mass_frac
    function = rho_fun
  []
[]

[AuxKernels]
  [rho]
    type = ADMaterialRealAux
    variable = rho
    property = rho
    execute_on = 'timestep_end'
  []
  [U_x]
    type = ADMaterialRealAux
    variable = U_x
    property = vel_x
    execute_on = 'timestep_end'
  []
  [U_y]
    type = ADMaterialRealAux
    variable = U_y
    property = vel_y
    execute_on = 'timestep_end'
  []
  [pressure]
    type = ADMaterialRealAux
    variable = pressure
    property = p
    execute_on = 'timestep_end'
  []
  [temperature]
    type = ADMaterialRealAux
    variable = temperature
    property = T_fluid
    execute_on = 'timestep_end'
  []
  [rho_et]
    type = ADMaterialRealAux
    variable = rho_et
    property = rho_et
    execute_on = 'timestep_end'
  []
  [courant]
    type = Courant
    variable = courant
    u = U_x
    v = U_y
  []
[]

[FVKernels]
  [mass_frac_time]
    type = PCNSFVDensityTimeDerivative
    variable = mass_frac
    rho = rho
  []
  [mass_frac_advection]
    type = PCNSFVKT
    variable = mass_frac
    eqn = "scalar"
  []
[]

[Functions]
  [ud_in]
    type = ParsedVectorFunction
    value_x = '0'
    value_y = '${v_in}'
  []
  [rho_fun]
    type = ParsedFunction
    value = '${rho_initial}*exp((x-0.2)^2)'
  []
[]

[FVBCs]
  [mf_bottom]
    type = PCNSFVStrongBC
    boundary = 'bottom'
    variable = mass_frac
    superficial_velocity = 'ud_in'
    T_fluid = ${T}
    eqn = 'scalar'
  []
[]

[Materials]
  [var_mat]
    type = PorousConservedVarMaterial
    rho = rho
    rho_et = rho_et
    superficial_rhou = rho_u
    superficial_rhov = rho_v
    fp = fp
    porosity = porosity
  []
  [porosity]
    type = GenericConstantMaterial
    prop_names = 'porosity'
    prop_values = '1'
  []
[]

[Executioner]
  solve_type = NEWTON
  type = Transient
  nl_max_its = 10
  steady_state_detection = true
  steady_state_tolerance = 1e-12
  abort_on_solve_fail = false
  dt = 5e-4
  num_steps = 100
[]

[Outputs]
  [out]
    type = Exodus
    execute_on = 'initial timestep_end'
  []
  [dof]
    type = DOFMap
    execute_on = 'initial'
  []
[]

[Debug]
  show_var_residual_norms = true
[]
