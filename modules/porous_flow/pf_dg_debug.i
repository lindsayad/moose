# transient
# using the PorousFlowFullySaturatedDarcyBase Kernel
[Mesh]
  [kgd_model_mesh]
    type = FileMeshGenerator
    file = two_element.e
  []
  [breakmesh]
    type = BreakMeshByBlockGenerator
    input = kgd_model_mesh
    block_pairs = 'X_Plus_Block X_Minus_Block'
    interface_name = DG_Interface
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [pp]
    initial_condition = 0.0
  []
[]

[Kernels]
  [mass0]
    #type = PorousFlowFullySaturatedMassTimeDerivative
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pp
  []
  [flux]
    type = PorousFlowFullySaturatedDarcyBase
    variable = pp
    gravity = '0 0 0'
  []
[]

[InterfaceKernels]
  [porepressure]
    type = PorousFlowFullySaturatedDarcyDG
    variable = pp
    neighbor_var = pp
    boundary = 'DG_Interface'
    stabilized_para = 1.0
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
[]

[Modules]
  [FluidProperties]
    [simple_fluid]
      type = SimpleFluidProperties
      bulk_modulus = 2e9
      density0 = 1000
      thermal_expansion = 0
      viscosity = 1e-3
    []
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
  []
  [ppss]
    type = PorousFlow1PhaseFullySaturated
    porepressure = pp
  []
  [massfrac]
    type = PorousFlowMassFraction
  []
  [simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.1
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1E-14 0 0 0 1E-14 0 0 0 1E-14'
  []
  [relperm]
    type = PorousFlowRelativePermeabilityCorey
    n = 0
    phase = 0
  []
  [solid_compilance]
    type = PorousFlowConstantBiotModulus
    solid_bulk_compliance = 1.46e-8
    biot_coefficient = 0.75
    block = 'X_Plus_Block X_Minus_Block'
  []
[]

[BCs]
  [Y_Plus]
    type = DirichletBC
    boundary = X_Plus
    value = 0
    variable = pp
  []
  [Injection]
    type = PorousFlowSink
    boundary = Center_Face
    variable = pp
    flux_function = -1 # a sink
  []
[]
[Preconditioning]
  #active = basic
  active = preferred
  [basic]
    type = SMP
    full = true
  []
  [preferred]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  []
[]

[Executioner]
  type = Transient
  solve_type = Newton
  l_max_its = 400
  l_tol = 1e-4
  nl_max_its = 200
  nl_abs_tol = 1e-6
  nl_rel_tol = 1e-4
  dt = 20
  end_time = 20
[]

[Outputs]
  file_base = debug_DG
  print_linear_residuals = false
  exodus = true
[]
