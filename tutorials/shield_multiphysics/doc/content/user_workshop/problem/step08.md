# Step 8: Mesh Adaptivity id=step08

!!end-intro

!---

## Step 8a: Coarse Solution

!listing step8a_coarse.i

!---

## Step 8a: Run

```bash
cd ~/projects/moose/tutorials/shield_multiphysics/step08_adaptivity
make -j 12 # use number of processors for your system
cd problems
../moose-opt -i step8a_coarse.i
```

!---

## Step 8b: Fine Solution

!listing step8b_fine.i

!---

## Step 8b: Run

```bash
cd ~/projects/moose/tutorials/shield_multiphysics/step08_adaptivity
make -j 12 # use number of processors for your system
cd problems
../moose-opt -i step8b_fine.i
```

!---

## Step 8c: Adaptive Mesh Solution

!listing step8c_adapt.i

!---

## Step 8c: Run

With the step 8 executable:

```bash
cd ~/projects/moose/tutorials/shield_multiphysics/step08_adaptivity
make -j 12 # use number of processors for your system
cd problems
../moose-opt -i step8c_adapt.i
```

With a conda MOOSE executable:

```bash
cd ~/projects/moose/tutorials/shield_multiphysics/step08_adaptivity/problems
conda activate moose
moose-opt -i step8c_adapt.i
```
