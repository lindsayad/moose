!include evbc_primary_element_ghosting.i

[Mesh]
  displacements = 'disp_x disp_y'
[]

[AuxVariables]
  [disp_x]
  []
  [disp_y]
  []
[]

[Constraints]
  [top]
    # Bind this constraint to the displaced mesh/problem so that the element it retains via
    # addRetainedGhostedElem() is kept alive only on the displaced DistributedMesh, while the
    # send-list augmentation that consumes the same ghosted-elements entry always resolves
    # against the undisplaced mesh.
    use_displaced_mesh = true
  []
[]

[Outputs]
  hide = 'disp_x disp_y'
[]
