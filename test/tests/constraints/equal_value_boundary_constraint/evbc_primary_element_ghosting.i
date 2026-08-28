[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 2
  []
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = u
  []
[]

[BCs]
  [neumann]
    type = NeumannBC
    variable = u
    boundary = 'left bottom'
    value = -1
  []
  [right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  []
[]

[UserObjects]
  # This UserObject introduces late geometric ghosting, which defers remote-element
  # deletion to the delete_remote_elements_after_late_geometric_ghosting task (run after
  # the constraint is added). That deletion is what removes the EVBC primary element unless
  # it is explicitly kept, so it is required to reproduce the bug.
  [late_geometric_ghosting]
    type = TestGhostBoundarySideUserObject
    boundary = left
  []
[]

[Constraints]
  [top]
    type = EqualValueBoundaryConstraint
    variable = u
    secondary = top
    primary_node_coord = '1 1 0'
    penalty = 1e7
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-ksp_norm_type -pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'preconditioned lu       mumps'
[]

[Outputs]
  exodus = true
[]
