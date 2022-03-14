mu = 1
rho = 1
k = 1
cp = 1
alpha = 1
# rayleigh=1e3
# hot_temp=${rayleigh}
hot_temp = 1
temp_ref=${fparse hot_temp / 2.}
l = 1

[Problem]
  kernel_coverage_check = false
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${l}
    ymin = 0
    ymax = ${l}
    nx = 8
    ny = 8
  []
  # second_order = true
[]

[Variables]
  [vel]
    family = LAGRANGE_VEC
    # order = SECOND
  []
  [pressure]
  []
  [T]
  []
  [lambda]
    family = SCALAR
    order = FIRST
  []
[]

[ICs]
  [vel]
    type = VectorConstantIC
    variable = vel
    x_value = 1e-12
    y_value = 1e-12
  []
[]

[Kernels]
  [mass]
    type = INSADMass
    variable = pressure
  []
  [mass_pspg]
    type = INSADMassPSPG
    variable = pressure
  []
  [mean_zero_pressure]
    type = ScalarLagrangeMultiplier
    variable = pressure
    lambda = lambda
  []

  [./momentum_viscous]
    type = INSADMomentumViscous
    variable = vel
  [../]
  [momentum_advection]
    type = INSADMomentumAdvection
    variable = vel
  []
  [momentum_pressure]
    type = INSADMomentumPressure
    variable = vel
    pressure = pressure
    integrate_p_by_parts = true
  []
  [./buoyancy]
    type = INSADBoussinesqBodyForce
    variable = vel
    temperature = T
    gravity = '0 -1 0'
  [../]
  # [./momentum_time]
  # type = INSADMomentumTimeDerivative
  # variable = vel
  # [../]
  [./gravity]
    type = INSADGravityForce
    variable = vel
    gravity = '0 -1 0'
  [../]
  [supg]
    type = INSADMomentumSUPG
    variable = vel
    velocity = vel
  []

  # [temperature_time]
  #   type = INSADHeatConductionTimeDerivative
  #   variable = T
  # []
  [temp_advection]
    type = INSADEnergyAdvection
    variable = T
  []
  [temp_conduction]
    type = ADHeatConduction
    variable = T
    thermal_conductivity = 'k'
  [../]
  [temp_supg]
    type = INSADEnergySUPG
    variable = T
    velocity = vel
  []
[]

[BCs]
  # Boundary conditions for laminar flow
  [wall]
    type = VectorDirichletBC
    variable = 'vel'
    boundary = 'top bottom left right'
    values = '0 0 0'
  []
  [T_hot]
    type = DirichletBC
    variable = T
    boundary = 'left'
    value = ${hot_temp}
  []
  [T_cold]
    type = DirichletBC
    variable = T
    boundary = 'right'
    value = 0
  []
[]

[Materials]
  [const_functor]
    type = ADGenericConstantMaterial
    prop_names = 'alpha cp k mu rho'
    prop_values = '${alpha} ${cp} ${k} ${mu} ${rho}'
  []
  [./const]
    type = GenericConstantMaterial
    prop_names =  'temp_ref'
    prop_values = '${temp_ref}'
  [../]
  [ins_mat]
    type = INSADStabilized3Eqn
    velocity = vel
    pressure = pressure
    temperature = T
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options = '-pc_svd_monitor'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'svd'
  # petsc_options_iname = '-pc_type -pc_factor_shift_type'
  # petsc_options_value = 'lu       NONZERO'

  nl_rel_tol = 1e-12
  # nl_abs_tol = 1e-8
  # nl_max_its = 10
  # line_search = 'none'
  # l_max_its = 20
  # end_time = 1e9
  # [TimeStepper]
  #   type = IterationAdaptiveDT
  #   optimal_iterations = 6
  #   dt = 1e-3
  #   growth_factor = 1.2
  # []
[]

[Outputs]
  exodus = true
  checkpoint = true
[]

[Postprocessors]
  [Rayleigh]
    type = RayleighNumber
    beta = ${alpha}
    T_hot = ${hot_temp}
    T_cold = 0
    rho_ave = ${rho}
    l = ${l}
    mu_ave = ${mu}
    k_ave = ${k}
    cp_ave = ${cp}
    gravity_magnitude = 1
  []
[]
