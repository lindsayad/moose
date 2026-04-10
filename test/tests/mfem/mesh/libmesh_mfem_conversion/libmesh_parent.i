[Mesh]
  type = FileMesh
  file = ../square.e
[]

[Problem]
  type = MFEMProblem
  first_order_mesh = false
[]

[FESpaces]
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
  []
[]

[Variables]
  [u]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[BCs]
  [sides]
    type = MFEMScalarDirichletBC
    variable = u
    coefficient = 1.0
  []
[]

[Kernels]
  [diff]
    type = MFEMDiffusionKernel
    variable = u
  []
[]

[Preconditioner]
  [boomeramg]
    type = MFEMHypreBoomerAMG
  []
[]

[Solver]
  type = MFEMHypreGMRES
  preconditioner = boomeramg
  l_tol = 1e-16
  l_max_its = 1000
[]

[Executioner]
  type = MFEMSteady
[]

[Outputs]
  file_base = 'libmesh_to_mfem_mesh_quads'
  csv = true
[]

[Postprocessors]
  [L2Error]
    type = MFEML2Error
    variable = u
    function = 1.0
  []
[]
