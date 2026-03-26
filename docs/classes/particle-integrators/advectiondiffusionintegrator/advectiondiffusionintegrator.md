---
layout: default
title: AdvectionDiffusionIntegrator
parent: AdvectionDiffusionIntegrator
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  AdvectionDiffusionIntegrator

Create an advection-diffusion integrator.


---

## Declaration
```matlab
 self = AdvectionDiffusionIntegrator(kinematicModel,kappa)
```
## Parameters
+ `kinematicModel`  `KinematicModel` instance defining `u(t,x,y)` and `v(t,x,y)`
+ `kappa`  scalar diffusivity in $$m^2 s^{-1}$$

## Returns
+ `self`  `AdvectionDiffusionIntegrator` instance

## Discussion

  `kinematicModel` supplies the deterministic velocity field,
  and `kappa` sets the scalar tracer diffusivity.


