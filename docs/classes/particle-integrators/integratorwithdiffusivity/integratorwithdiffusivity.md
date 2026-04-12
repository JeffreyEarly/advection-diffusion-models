---
layout: default
title: IntegratorWithDiffusivity
parent: IntegratorWithDiffusivity
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  IntegratorWithDiffusivity

Create an additive-diffusion integrator on a box domain.


---

## Declaration
```matlab
 self = IntegratorWithDiffusivity(f,y0,dt=...,kappa=...,ymin=...,ymax=...)
```
## Parameters
+ `f`  deterministic drift function $$f(t,y)$$
+ `y0`  initial condition $$y_0$$ stored as `nParticles x nDims`
+ `dt`  positive scalar timestep $$dt$$
+ `kappa`  scalar or `1 x nDims` diffusivity vector $$\kappa$$; defaults to `0`
+ `ymin`  scalar or `1 x nDims` lower box bounds; defaults to `-Inf`
+ `ymax`  scalar or `1 x nDims` upper box bounds; defaults to `Inf`

## Discussion

  The constructor keeps the legacy scalar-expansion behavior
  for `kappa`, `ymin`, and `ymax`, with defaults corresponding
  to zero diffusivity on an unbounded domain.
