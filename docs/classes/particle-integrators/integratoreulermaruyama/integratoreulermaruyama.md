---
layout: default
title: IntegratorEulerMaruyama
parent: IntegratorEulerMaruyama
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  IntegratorEulerMaruyama

Create an Euler-Maruyama integrator for $$dy = f\,dt + g\,dW_t$$.


---

## Declaration
```matlab
 self = IntegratorEulerMaruyama(f,g,y0,dt=...)
```
## Parameters
+ `f`  deterministic drift function $$f(t,y)$$
+ `g`  stochastic amplitude function $$g(t,y)$$
+ `y0`  initial condition $$y_0$$ stored as `nParticles x nDims`
+ `dt`  positive scalar timestep $$dt$$

## Discussion

  The constructor validates both `f(0,y0)` and `g(0,y0)` and
  requires them to return arrays with the same shape as `y0`.


