---
layout: default
title: IntegratorWithObstacles
parent: IntegratorWithObstacles
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  IntegratorWithObstacles

Create a 2D obstacle-aware stochastic integrator.


---

## Declaration
```matlab
 self = IntegratorWithObstacles(f,y0,dt=...,kappa=...,ymin=...,ymax=...,obstacles=polyshape.empty(0,1),isPeriodic=false)
```
## Parameters
+ `f`  deterministic drift function $$f(t,y)$$
+ `y0`  initial condition $$y_0$$ stored as `nParticles x 2`
+ `dt`  positive scalar timestep $$dt$$
+ `kappa`  scalar or `1 x 2` diffusivity vector $$\kappa$$; defaults to `0`
+ `ymin`  scalar or `1 x 2` lower coordinate bounds; defaults to `-Inf`
+ `ymax`  scalar or `1 x 2` upper coordinate bounds; defaults to `Inf`
+ `obstacles`  non-overlapping reflecting `polyshape` obstacles
+ `isPeriodic`  logical scalar or `1 x 2` periodicity flags; defaults to `false`

## Discussion

  The constructor validates the 2D state layout, expands
  scalar-valued `kappa`, `ymin`, `ymax`, and `isPeriodic`,
  stores the domain bounds, and precomputes polygon and
  bounding-box geometry used by the reflection algorithm.


