<!-- MOOSE Documentation Stub: Remove this when content is added. -->

# UserForcingFunction

## Description

`UserForcingFunction` is a more restricted version of
[`BodyForce`](systems/Kernels/framework/BodyForce.md) that implements a body force or
source term as a function of space and/or time. The forcing function is provided
through the parameter `function`.

## Example Syntax

The block below shows `UserForcingFunction` used in the context of a
transient-diffusion-reaction problem. The `f_fn` object is assigned to the
`function` parameter.

!listing test/tests/kernels/ode/ode_sys_impl_test.i block=Kernels label=false

!syntax parameters /Kernels/UserForcingFunction

!syntax inputs /Kernels/UserForcingFunction

!syntax children /Kernels/UserForcingFunction
