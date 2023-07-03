[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 2
  ny = 2
[]

[Variables]
  [u]
  []
  [v]
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = u
  []
  [diff_v]
    type = Diffusion
    variable = v
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  []
  [right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  []
  [left_v]
    type = DirichletBC
    variable = v
    boundary = left
    value = 0
  []
  [right_v]
    type = DirichletBC
    variable = v
    boundary = right
    value = 1
  []
[]

[Preconditioning]
  [FSP]
    type = FSP
    topsplit = 'by_var'
    [by_var]
      splitting = 'u v'
      splitting_type  = additive
      petsc_options_iname = '-ksp_type'
      petsc_options_value = 'fgmres'
    []
      [u]
        vars = 'u'
        splitting = 'u_diri u_bulk'
        splitting_type = multiplicative
        petsc_options_iname = '-ksp_type'
        petsc_options_value = 'fgmres'
      []
        [u_diri]
          vars = 'u'
          petsc_options = '-ksp_view_pmat'
          sides = 'left right'
        []
        [u_bulk]
          vars = 'u'
          petsc_options = '-ksp_view_pmat'
          petsc_options_iname = '-ksp_type'
          petsc_options_value = 'cg'
          unsides = 'left right'
        []
      [v]
        vars = 'v'
        splitting = 'v_diri v_bulk'
        splitting_type = multiplicative
        petsc_options_iname = '-ksp_type'
        petsc_options_value = 'fgmres'
      []
        [v_diri]
          vars = 'v'
          petsc_options = '-ksp_view_pmat'
          sides = 'left right'
        []
        [v_bulk]
          vars = 'v'
          petsc_options = '-ksp_view_pmat'
          petsc_options_iname = '-ksp_type'
          petsc_options_value = 'cg'
          unsides = 'left right'
        []
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
[]

[Outputs]
  exodus = true
[]
