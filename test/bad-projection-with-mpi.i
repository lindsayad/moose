[Mesh]
        [fmg]
                type = FileMeshGenerator
                file = 'refined_mesh.e' 
                skip_partitioning = False
                use_for_exodus_restart = False
        []
[]
[AuxVariables]
        [./TK_elem]
                order = CONSTANT
                family = MONOMIAL
        [../]
        [./TK_proj]
                order = FIRST
                family = LAGRANGE
        [../]
[]
[AuxKernels]
        [./update_TK_elem]
            type = FunctionAux
            variable = TK_elem
            function = aux_TK
            execute_on = 'INITIAL TIMESTEP_END'
        [../]
        [./update_TK_proj]
            type = ProjectionAux
            v = TK_elem
            execute_on = 'TIMESTEP_END'
            variable = TK_proj
        [../]
[]
[Functions]
        [./aux_TK]
                type = ParsedFunction
                expression = 'Tmin + (Tmax - Tmin)*(x+y+z-xmin-ymin-zmin)/(xmax-xmin+ymax-ymin+zmax-zmin)'
                symbol_names = 'Tmin   Tmax     xmin  xmax     ymax   ymin  zmax   zmin'
                symbol_values= '300.0  2000.0   0.0   2.8e-3   8e-4   0.0   2e-3   0.0'
        [../]
[]
[Outputs]
   [./out]
           type = Exodus
           execute_on = 'final'
   [../]
[]
[Executioner]
  type = Transient
   dt = 1.0
  end_time = 2.0
[]
[Variables]
        [./u]
        [../]
[]
[Kernels]
        [./null]
                type = NullKernel
                variable = u
        [../]
[]