mu = 1e-2
rho = 1
l = 2
U = 1
advected_interp_method = 'upwind'
velocity_interp_method = 'rc'

[GlobalParams]
  rhie_chow_user_object = 'rc'
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = vel_x
    v = vel_y
    pressure = pressure
  []
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 10
    ymin = ${fparse -l / 2}
    ymax = ${fparse l / 2}
    nx = 100
    ny = 20
  []
  uniform_refine = 0
[]

[Variables]
  [vel_x]
    type = INSFVVelocityVariable
    initial_condition = 0
  []
  [vel_y]
    type = INSFVVelocityVariable
    initial_condition = 0
  []
  [pressure]
    type = INSFVPressureVariable
  []
[]

[FVKernels]
  [mass]
    type = INSFVMassAdvection
    variable = pressure
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
  []

  [u_time]
    type = INSFVMomentumTimeDerivative
    variable = vel_x
    rho = ${rho}
    momentum_component = 'x'
  []
  [u_advection]
    type = INSFVMomentumAdvection
    variable = vel_x
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
    momentum_component = 'x'
  []
  [u_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_x
    mu = ${mu}
    momentum_component = 'x'
  []
  [u_pressure]
    type = INSFVMomentumPressure
    variable = vel_x
    momentum_component = 'x'
    pressure = pressure
  []

  [v_time]
    type = INSFVMomentumTimeDerivative
    variable = vel_y
    rho = ${rho}
    momentum_component = 'y'
  []
  [v_advection]
    type = INSFVMomentumAdvection
    variable = vel_y
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
    momentum_component = 'y'
  []
  [v_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_y
    mu = ${mu}
    momentum_component = 'y'
  []
  [v_pressure]
    type = INSFVMomentumPressure
    variable = vel_y
    momentum_component = 'y'
    pressure = pressure
  []
[]

[FVBCs]
  [inlet-u]
    type = INSFVInletVelocityBC
    boundary = 'left'
    variable = vel_x
    function = '${U}'
  []
  [inlet-v]
    type = INSFVInletVelocityBC
    boundary = 'left'
    variable = vel_y
    function = '0'
  []
  [walls-u]
    type = INSFVNoSlipWallBC
    boundary = 'top bottom'
    variable = vel_x
    function = 0
  []
  [walls-v]
    type = INSFVNoSlipWallBC
    boundary = 'top bottom'
    variable = vel_y
    function = 0
  []
  [outlet_p]
    type = INSFVOutletPressureBC
    boundary = 'right'
    variable = pressure
    function = '0'
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  nl_abs_tol = 1e-12
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    dt = 1e-3
  []
  steady_state_detection = true
[]

[Preconditioning]
  [FSP]
    type = FSP
    topsplit = 'up'
    [up]
      splitting = 'u p'
      splitting_type  = schur
      petsc_options_iname = '-pc_fieldsplit_schur_fact_type  -pc_fieldsplit_schur_precondition -ksp_gmres_restart -ksp_rtol -ksp_type'
      petsc_options_value = 'full                            self                             300                1e-5      fgmres'
    []
    [u]
      vars = 'vel_x vel_y'
      petsc_options_iname = '-pc_type -pc_hypre_type -ksp_pc_side -ksp_type -ksp_rtol -pc_hypre_boomeramg_relax_type_all -pc_hypre_boomeramg_restriction_type -pc_hypre_boomeramg_postrelax -pc_hypre_boomeramg_grid_sweeps_up -pc_hypre_boomeramg_grid_sweeps_down'
      petsc_options_value = 'hypre    boomeramg      right        gmres     1e-2      SOR/Jacobi                         2                                    F,F,C                         3 0'
    []
    [p]
      vars = 'pressure'
      petsc_options = '-pc_lsc_scale_diag -ksp_monitor'
      petsc_options_iname = '-ksp_type -ksp_gmres_restart -ksp_rtol -pc_type -ksp_pc_side -pc_type  -lsc_pc_type -lsc_pc_factor_mat_solver_type -lsc_ksp_type -lsc_ksp_rtol'
      petsc_options_value = 'fgmres    300                1e-3      lsc      right        lsc       lu           mumps                           gmres        1e-1'
    []
  []
[]

[Outputs]
  print_linear_residuals = true
  print_nonlinear_residuals = true
  [out]
    type = Exodus
    hide = 'Re lin cum_lin'
  []
[]

[Postprocessors]
  [Re]
    type = ParsedPostprocessor
    function = '${rho} * ${l} * ${U}'
    pp_names = ''
  []
  [lin]
    type = NumLinearIterations
  []
  [cum_lin]
    type = CumulativeValuePostprocessor
    postprocessor = lin
  []
[]
