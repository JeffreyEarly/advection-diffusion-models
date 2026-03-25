---
layout: default
title: Integrator
parent: Integrator
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  Integrator

Create an RK4 integrator for $$\frac{dy}{dt} = f(t,y)$$.


---

## Declaration
```matlab
 self = Integrator(f,y0,dt=...)
```
## Parameters
+ `f`  function handle for the deterministic right-hand side $$f(t,y)$$
+ `y0`  initial condition $$y_0$$ stored as `nParticles x nDims`
+ `dt`  positive scalar timestep $$dt$$

## Discussion

  The constructor validates the initial right-hand side
  evaluation `f(0,y0)` and stores the fixed timestep `dt`. The
  state array `y0` is used without reshaping, so its
  `nParticles x nDims` layout becomes the integration contract
  for all later calls to `f`.


