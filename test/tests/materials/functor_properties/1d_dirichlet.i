[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 10
  xmax = 2
[]

[GlobalParams]
  mat_prop_name = 'test'
[]

[Variables]
  [v]
    family = MONOMIAL
    order = CONSTANT
    fv = true
  []
[]

[AuxVariables]
  [sink]
    type = MooseVariableFVReal
  []
[]

[ICs]
  [sink]
    type = FunctionIC
    variable = sink
    function = 'x^3'
  []
[]


[FVKernels]
  [diff]
    type = FVDiffusion
    variable = v
    coeff = 1
  []
  [sink]
    type = FVFunctorElementalKernel
    variable = v
  []
[]

[FVBCs]
  [bounds]
    type = FVDirichletBC
    variable = v
    boundary = 'left right'
    value = 0
  []
[]

[Materials]
  [functor]
    type = FVVarFunctorMaterial
    var = sink
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
[]

[Outputs]
  exodus = true
[]