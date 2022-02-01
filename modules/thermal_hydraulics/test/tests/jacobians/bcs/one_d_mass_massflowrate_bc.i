[JacobianTest1Phase]
  A = 1e-1
  p = 1e6
  T = 300
  vel = 2
  snes_test_err = 1e-4
  fp_1phase = fp_1phase
[]

[Modules/FluidProperties]
  [fp_1phase]
    type = StiffenedGasFluidProperties
    gamma = 2.35
    q = -1167e3
    q_prime = 0
    p_inf = 1.e9
    cv = 1816
  []
[]

[BCs]
  [bc_1]
    type = OneDMassMassFlowRateTemperatureBC
    variable = rhouA
    boundary = 0
    normal = -1
    m_dot = 1234
  []
  [bc_2]
    type = OneDMassMassFlowRateTemperatureBC
    variable = rhouA
    boundary = 1
    normal = 1
    m_dot = 1234
  []
[]