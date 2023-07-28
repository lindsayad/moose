mu = 1
rho = 1
n = 64

[GlobalParams]
  velocity_interp_method = 'average'
  advected_interp_method = 'upwind'
  rhie_chow_user_object = 'rc'
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    nx = ${n}
    ny = ${n}
  []
[]

[Variables]
  [vel_x]
    type = INSFVVelocityVariable
  []
  [vel_y]
    type = INSFVVelocityVariable
  []
  [pressure]
    type = INSFVPressureVariable
  []
[]

[AuxVariables]
  [U]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
[]

[AuxKernels]
  [mag]
    type = VectorMagnitudeAux
    variable = U
    x = vel_x
    y = vel_y
  []
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = vel_x
    v = vel_y
    pressure = pressure
  []
[]

[FVKernels]
  [mass_advection]
    type = INSFVMassAdvection
    variable = pressure
    rho = ${rho}
  []
  # [u_time]
  #   type = INSFVMomentumTimeDerivative
  #   variable = vel_x
  #   rho = ${rho}
  #   momentum_component = 'x'
  # []
  [u_advection]
    type = INSFVMomentumAdvection
    variable = vel_x
    rho = ${rho}
    momentum_component = 'x'
  []
  [u_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_x
    mu = 'mu'
    momentum_component = 'x'
  []
  [u_pressure]
    type = INSFVMomentumPressure
    variable = vel_x
    momentum_component = 'x'
    pressure = pressure
  []
  # [u_friction]
  #   type = INSFVMomentumFriction
  #   variable = vel_x
  #   momentum_component = 'x'
  #   linear_coef_name = 1
  # []

  # [v_time]
  #   type = INSFVMomentumTimeDerivative
  #   variable = vel_y
  #   rho = ${rho}
  #   momentum_component = 'y'
  # []
  [v_advection]
    type = INSFVMomentumAdvection
    variable = vel_y
    rho = ${rho}
    momentum_component = 'y'
  []
  [v_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_y
    mu = 'mu'
    momentum_component = 'y'
  []
  [v_pressure]
    type = INSFVMomentumPressure
    variable = vel_y
    momentum_component = 'y'
    pressure = pressure
  []
  # [v_friction]
  #   type = INSFVMomentumFriction
  #   variable = vel_y
  #   momentum_component = 'y'
  #   linear_coef_name = 1
  # []
[]

[FVBCs]
  [top_x]
    type = INSFVNoSlipWallBC
    variable = vel_x
    boundary = 'top'
    function = 1
  []
  [no_slip_x]
    type = INSFVNoSlipWallBC
    variable = vel_x
    boundary = 'left right bottom'
    function = 0
  []
  [no_slip_y]
    type = INSFVNoSlipWallBC
    variable = vel_y
    boundary = 'left right top bottom'
    function = 0
  []
[]

[Materials]
  [mu]
    type = ADGenericFunctorMaterial
    prop_names = 'mu'
    prop_values = '${mu}'
  []
[]

[Preconditioning]
  [FSP]
    type = FSP
    topsplit = 'up'
    [up]
      splitting = 'u p'
      splitting_type  = schur
      petsc_options_iname = '-pc_fieldsplit_schur_fact_type  -pc_fieldsplit_schur_precondition -ksp_gmres_restart -ksp_rtol -ksp_type -ksp_atol'
      petsc_options_value = 'full                            self                             300                1e-5      fgmres 1e-9'
    []
      [u]
        vars = 'vel_x vel_y'
        petsc_options_iname = '-pc_type'
        petsc_options_value = 'lu'
      []
      [p]
        vars = 'pressure'
        petsc_options = '-ksp_monitor -pc_lsc_scale_diag'
        petsc_options_iname = '-ksp_type -ksp_gmres_restart -ksp_rtol -pc_type -ksp_pc_side -lsc_pc_type -lsc_ksp_rtol -lsc_ksp_type'
        petsc_options_value = 'fgmres     300                1e-2     lsc      right        hypre        1e-1 gmres'
      []
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  # [TimeStepper]
  #   type = IterationAdaptiveDT
  #   optimal_iterations = 6
  #   dt = 1e-1
  # []
  # steady_state_detection = true
[]

[Outputs]
  exodus = true
[]
