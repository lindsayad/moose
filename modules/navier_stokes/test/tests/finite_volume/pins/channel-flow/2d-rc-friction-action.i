mu = 1.1
rho = 1

[Mesh]
  [mesh]
    type = CartesianMeshGenerator
    dim = 2
    dx = '2.5 2.5'
    dy = '1.0'
    ix = '20 20'
    iy = '20'
    subdomain_id = '1 2'
  []
[]

[Modules]
  [NavierStokesFV]
    compressibility = 'incompressible'
    porous_medium_treatment = true

    density = 'rho'
    dynamic_viscosity = 'mu'
    porosity = 'porosity'

    initial_velocity = '1 1e-6 0'
    initial_pressure = 0.0

    inlet_boundaries = 'left'
    momentum_inlet_types = 'fixed-velocity'
    momentum_inlet_function = '1 0'
    wall_boundaries = 'top bottom'
    momentum_wall_types = 'noslip symmetry'
    outlet_boundaries = 'right'
    momentum_outlet_types = 'fixed-pressure'
    pressure_function = '0'

    friction_types = 'darcy forchheimer'
    friction_coeffs = 'Darcy_coefficient Forchheimer_coefficient'

    mass_advection_interpolation = 'average'
    momentum_advection_interpolation = 'average'
  []
[]

[Materials]
  [const]
    type = ADGenericFunctorMaterial
    prop_names = 'rho mu'
    prop_values = '${rho} ${mu}'
  []
  [friction_factors]
    type = ADGenericVectorFunctorMaterial
    prop_names = 'Darcy_coefficient Forchheimer_coefficient'
    prop_values = 'D_old_convention D_old_convention D_old_convention
                   F_old_convention F_old_convention F_old_convention'
  []
  # These two materials are only used to adapt to the old convention for friction coefficient
  # This should never be used outside of this input file.
  [legacy_darcy_conversion]
    type = ADParsedFunctorMaterial
    property_name = 'D_old_convention'
    expression = '0.1 * rho / porosity / mu'
    functor_symbols = 'rho porosity mu'
    functor_names = '${rho} porosity ${mu}'
  []
  [legacy_forchheimer_conversion]
    type = ADParsedFunctorMaterial
    property_name = 'F_old_convention'
    expression = '0.2 / porosity / speed'
    functor_symbols = 'speed porosity'
    functor_names = 'speed porosity'
  []
[]

[AuxVariables]
  [porosity]
    family = MONOMIAL
    order = CONSTANT
    fv = true
    initial_condition = 0.5
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
  nl_rel_tol = 1e-11
  nl_abs_tol = 1e-14
[]

# Some basic Postprocessors to visually examine the solution
[Postprocessors]
  [inlet-p]
    type = SideIntegralVariablePostprocessor
    variable = pressure
    boundary = 'left'
  []
  [outlet-u]
    type = SideIntegralVariablePostprocessor
    variable = superficial_vel_x
    boundary = 'right'
  []
[]

[Outputs]
  exodus = true
[]
