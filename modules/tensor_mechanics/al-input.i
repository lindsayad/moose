[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = true
  # scaling = 1e0
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 5
  ny = 1
  nz = 1
  xmin = 0.0
  xmax = 4.0
  ymin = 0.0
  ymax = 0.5
  zmin = 0.0
  zmax = 1.0
  elem_type = HEX8
[]

[Variables]
  [disp_x]
    # scaling = 3.25447e-06
  []
  [disp_y]
    # scaling = 2.29038e-06
  []
  [disp_z]
    # scaling = 3.60462e-06
  []
[]

[Modules/TensorMechanics/Master]
  [./all]
    strain = FINITE
    generate_output = 'stress_xx stress_xy stress_yy strain_xx strain_xy strain_yy'
  [../]
[]

[Functions]
  [./loading_func]
    type = PiecewiseLinear
    x = '0.  5.'
    y = '0. 50.0'
  [../]
[]

[BCs]
  [./free_end_moment]
    type = Pressure
    variable = disp_y
    component = 0
    boundary = right
    factor = 1
    function = loading_func
  [../]
  [./FixedCenterLineX]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./FixedCenterLineY]
    type = DirichletBC
    variable = disp_y
    boundary = left
    value = 0.0
  [../]
  [./FixedCenterLineZ]
    type = DirichletBC
    variable = disp_z
    boundary = left
    value = 0.0
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e6
    poissons_ratio = 0.3
  [../]
  [./stress]
    type = ComputeFiniteStrainElasticStress
  [../]
[]

[Preconditioning]
  [./pc]
    type = SMP
    full = True
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu     superlu_dist'

  l_max_its = 20
  l_tol = 1e-2

  nl_max_its = 50
  # nl_rel_tol = 1e-8
  # nl_abs_tol = 1e-4

  start_time = 0.0
  dt = 1
  end_time = 5
[]

[Postprocessors]
  [./num_lin_it]
    type = NumLinearIterations
  [../]
  [./num_nonlin_it]
    type = NumNonlinearIterations
  [../]
  [./tot_lin_it]
    type = CumulativeValuePostprocessor
    postprocessor = num_lin_it
  [../]
  [./tot_nonlin_it]
    type = CumulativeValuePostprocessor
    postprocessor = num_nonlin_it
  [../]
  [./alive_time]
    type = PerfGraphData
    section_name = Root
    data_type = TOTAL
  [../]
  [./max_beam_deflection]
    type = NodalMaxValue
    variable = disp_y
    boundary = 'right'
  [../]
[]

[Outputs]
  exodus = true
  csv = true
[]
