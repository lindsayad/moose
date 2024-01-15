# Step 9: Postprocessors id=step09

!---

Aggregate values based on simulation data are useful for understanding the simulation as well
as defining coupling values across coupled equations.

There are two main systems for aggregating data: Postprocessors and VectorPostprocessors.

!!end-intro

!---

## Step 9: Input File

!listing step09_postprocessors/inputs/step9.i

!---

## Step 9: Run

```bash
cd ~/projects/moose/tutorials/shield_multiphysics/step09_postprocessors
make -j 12 # use number of processors for your system
cd problems
../moose-opt -i step9.i
```
