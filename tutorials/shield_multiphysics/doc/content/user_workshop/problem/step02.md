# Step 2: Simple Heat Conduction Kernel id=step02

!---

## Kernel Object

To implement the Heat Conduction equation, a `Kernel` object is needed to add the conductivity coefficient
to the diffusion equation.

!equation
-\nabla \cdot k \nabla T = 0,

where $k$ is the thermal diffusivity.

A `Kernel` is C++ class, which inherits from `MooseObject` that is used by MOOSE for coding
volume integrals of a [!ac](PDE).

!!end-intro

!---

## CoefDiffusion Kernel

To implement the coefficient a new Kernel object must be created: `CoefDiffusion`.

This object will inherit from Diffusion and will use input parameters for specifying the
diffusivity.

!---

## CoefDiffusion.h

!listing step02_coef_diffusion/include/kernels/CoefDiffusion.h

!---

## CoefDiffusion.C

!listing step02_coef_diffusion/src/kernels/CoefDiffusion.C

!---

## Step 2: Input File

We introduce block restriction to differentiate between water and concrete.

!listing step02_coef_diffusion/problems/step2.i

!---

## Step 2: Run

```bash
cd ~/projects/moose/tutorials/shield_multiphysics/step02_coef_diffusion
make -j 12 # use number of processors for your system
cd problems
../moose-opt -i step2.i
```

!---

## Step 2: Result

!media shield_multiphysics/results/step02.png style=width:70%;
