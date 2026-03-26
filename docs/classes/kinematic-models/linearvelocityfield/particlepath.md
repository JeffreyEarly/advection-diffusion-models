---
layout: default
title: particlePath
parent: LinearVelocityField
grand_parent: Classes
nav_order: 4
mathjax: true
---

#  particlePath

Evaluate the analytical particle trajectory solution.


---

## Declaration
```matlab
 [x, y] = particlePath(self,x0,y0,t,kappa,u_0,v_0)
```
## Parameters
+ `x0`  initial x positions in meters
+ `y0`  initial y positions in meters
+ `t`  output time vector in seconds
+ `kappa`  scalar diffusivity in $$m^2 s^{-1}$$
+ `u_0`  background x-velocity in $$m s^{-1}$$ used by the analytical formulas
+ `v_0`  background y-velocity in $$m s^{-1}$$ used by the analytical formulas

## Returns
+ `x`  x positions in meters with shape `[length(t) nParticles]`
+ `y`  y positions in meters with shape `[length(t) nParticles]`

## Discussion

  Initial particle positions are normalized to row vectors, and the output
  time vector is normalized to a column vector before the analytical
  update is applied.


