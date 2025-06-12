#Toy Geometry 2 to test mortar contact.

[Mesh]

  type = MeshGeneratorMesh
  patch_update_strategy = iteration

  #
  #Binder geometry
  #

  [binder_block_tg23]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    ny = 10
    xmin = 0
    xmax = 2e+2
    ymin = 0
    ymax = 2e+2
    boundary_id_offset = 800
    subdomain_ids = 800
    boundary_name_prefix = 'binder'
    elem_type = QUAD4
  []

  [binder_rm]
    type = SubdomainBoundingBoxGenerator
    input = binder_block_tg23
    block_id = 200
    bottom_left = '1e+2 0 0'
    top_right = '2e+2 1e+2 0'
  []

  [binder_rm_nodeset]
    type = BoundingBoxNodeSetGenerator
    input = binder_rm
    new_boundary = 801
    bottom_left = '1e+2 0 0'
    top_right = '2e+2 1e+2 0'
  []

  [binder]
    type = BlockDeletionGenerator
    block = 200
    input = binder_rm_nodeset
    new_boundary = 815
  []

  [binder_inner_top]
    type = SideSetsFromBoundingBoxGenerator
    input = binder
    included_boundaries = 815
    boundary_new = 812
    bottom_left = '0.95e+2 0.9e+2 0'
    top_right = '2.1e+2 1.1e+2 0'
  []

  [binder_inner_left]
    type = SideSetsFromBoundingBoxGenerator
    input = binder
    included_boundaries = 815
    boundary_new = 813
    bottom_left = '0.9e+2 0 0'
    top_right = '1.1e+2 1.05e+2 0'
  []

  #
  #Crystal geometry
  #

  [crystal]
    type = GeneratedMeshGenerator
    nx = 10
    ny = 9
    dim = 2
    xmin = 1e+2
    xmax = 2e+2
    ymin = 0
    ymax = 1e+2
    boundary_id_offset = 100
    subdomain_ids = 100
    boundary_name_prefix = 'sd1'
  []

  [crystal_subdomain]
    type = MeshCollectionGenerator
    inputs = 'binder binder_inner_top binder_inner_left crystal' # bcit_boundary ct_boundary'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [disp_x]
    order = FIRST
    family = LAGRANGE
  []
  [disp_y]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]

  #Initialize independent variables

  [temp]
    order = CONSTANT
    family = MONOMIAL
  []
  [irr]
    order = CONSTANT
    family = MONOMIAL
  []
  [initial_x]
    order = CONSTANT
    family = MONOMIAL
  []

  #Initialize dependent variables

  [strain_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_y]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Functions]

  [temp_def]
    type = ConstantFunction
    value = 800
  []
  [irr_def]
    type = ConstantFunction
    value = 2
  []
  [initial_x_def]
    type = ConstantFunction
    value = 100
  []
[]

[Physics]

  [SolidMechanics]

    [QuasiStatic]
      [all]
        eigenstrain_names = 'thermal_strain irr_strain'
        add_variables = true
        generate_output = 'vonmises_stress'
        block = '100 800'
      []
    []
  []
[]

[Contact]

  [bc_inner_top]
    primary = 812
    secondary = 102
    formulation = mortar
    model = frictionless
  []
[]

[AuxKernels]

  [initial_x]
    type = FunctionAux
    variable = initial_x
    function = initial_x_def
    use_displaced_mesh = false
  []
  [tempfuncaux]
    type = FunctionAux
    variable = temp
    function = temp_def
    use_displaced_mesh = false
  []
  [irrfuncaux]
    type = FunctionAux
    variable = irr
    function = irr_def
    use_displaced_mesh = false
  []
[]

[BCs]

  [right]
    type = DirichletBC
    boundary = 'binder_right sd1_right'
    variable = disp_x
    value = 0.
  []
  [bottom]
    type = DirichletBC
    boundary = 'binder_bottom sd1_bottom'
    variable = disp_y
    value = 0.
  []
[]

[Materials]

  #
  # Binder properties - isotropic elasticity + 2 eigenstrains
  #

  [binder_elasticity_tensor]
    block = '800 801'
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 10e6
    poissons_ratio = 0.01
  []
  [binder_therm_prefactor]
    type = DerivativeParsedMaterial
    block = '800 801'
    coupled_variables = 'temp'
    property_name = binder_therm_prefactor
    constant_names = 'a T'
    constant_expressions = '1.3e-5 298'
    expression = '(a*(temp-T))'
  []
  [binder_thermal_strain]
    type = ComputeVariableEigenstrain
    block = '800 801'
    eigen_base = '1 0 0 0 1 0 0 0 1'
    args = 'temp'
    prefactor = binder_therm_prefactor
    eigenstrain_name = thermal_strain
  []
  [binder_irr_prefactor]
    type = DerivativeParsedMaterial
    block = '800 801'
    coupled_variables = 'irr initial_x'
    property_name = binder_irr_prefactor
    constant_names = 'm'
    constant_expressions = '0'
    expression = '((m*irr)/100)'
  []
  [binder_irr_strain]
    type = ComputeVariableEigenstrain
    block = '800 801'
    eigen_base = '1 0 0 0 1 0 0 0 1'
    args = 'irr'
    prefactor = binder_irr_prefactor
    eigenstrain_name = irr_strain
  []

  #
  # Crystal properties - orthotropic elasticity + 2 eigenstrains
  #

  [elasticity_tensor]
    type = ComputeElasticityTensor
    block = '100'
    fill_method = orthotropic
    C_ijkl = '1.095e12 3.65e10 1.095e12 2.8568e8 9.549e6 9.549e6 0.01 0.01 0.3 0.3 0.01 0.01'
  []
  [therm_prefactor]
    type = DerivativeParsedMaterial
    block = '100'
    coupled_variables = 'temp'
    property_name = therm_prefactor
    constant_names = 'a T'
    constant_expressions = '1.3e-5 298'
    expression = '(a*(temp-T))'
  []
  [thermal_strain]
    type = ComputeVariableEigenstrain
    block = '100'
    #eigen_base = '1 0 0 0 1 0 0 0 1'
    eigen_base = '-0.0577 0 0 0 1 0 0 0 1'
    args = 'temp'
    prefactor = therm_prefactor
    eigenstrain_name = thermal_strain
  []
  [irr_prefactor]
    type = DerivativeParsedMaterial
    block = '100'
    coupled_variables = 'irr initial_x'
    property_name = irr_prefactor
    constant_names = 'm'
    constant_expressions = '1.185'
    expression = '((m*irr)/100)'
  []
  [irr_strain]
    type = ComputeVariableEigenstrain
    block = '100'
    #eigen_base = '1 0 0 0 1 0 0 0 1'
    eigen_base = '-0.31 0 0 0 1 0 0 0 1'
    args = 'irr'
    prefactor = irr_prefactor
    eigenstrain_name = irr_strain
  []

  [stress]
    type = ComputeLinearElasticStress
    block = '100 800'
  []
[]

[Preconditioning]
  [prec1]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Steady
  solve_type = 'Newton'
  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'
  nl_rel_tol = 0
  nl_abs_tol = 1e-2
[]

[Outputs]
  exodus = true
[]
