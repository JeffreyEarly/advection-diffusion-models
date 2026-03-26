---
layout: default
title: particleTrajectoriesWithStochasticVelocity
parent: AdvectionDiffusionIntegrator
grand_parent: Classes
nav_order: 5
mathjax: true
---

#  particleTrajectoriesWithStochasticVelocity

Integrate particles with extra stochastic velocity terms.


---

## Declaration
```matlab
 [t, x, y] = particleTrajectoriesWithStochasticVelocity(self,x0,y0,T,dt,u,v)
```
## Parameters
+ `x0`  initial x positions in meters; any input shape is flattened to a column vector
+ `y0`  initial y positions in meters; any input shape is flattened to a column vector
+ `T`  total integration duration in seconds
+ `dt`  output time increment in seconds
+ `u`  function handle `u(t)` returning one stochastic x-velocity value per particle
+ `v`  function handle `v(t)` returning one stochastic y-velocity value per particle

## Returns
+ `t`  column vector of output times in seconds
+ `x`  x positions in meters with shape `[length(t) nParticles]`
+ `y`  y positions in meters with shape `[length(t) nParticles]`

## Discussion

  This method integrates

  $$ d\mathbf{r}(t) = [\mathbf{u}_{model}(t,\mathbf{r}) + \mathbf{u}_{stochastic}(t)]\,dt + \sqrt{2\kappa}\,d\mathbf{W}_t, $$

  where `u(t)` and `v(t)` return column vectors with one value per active
  particle at the current output time.


