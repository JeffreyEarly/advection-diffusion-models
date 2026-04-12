---
layout: default
title: particleTrajectories
parent: AdvectionDiffusionIntegrator
grand_parent: Classes
nav_order: 4
mathjax: true
---

#  particleTrajectories

Integrate particles from initial positions.


---

## Declaration
```matlab
 [t, x, y] = particleTrajectories(self,x0,y0,T,dt)
```
## Parameters
+ `x0`  initial x positions in meters; any input shape is flattened to a column vector
+ `y0`  initial y positions in meters; any input shape is flattened to a column vector
+ `T`  total integration duration in seconds
+ `dt`  output time increment in seconds

## Returns
+ `t`  column vector of output times in seconds
+ `x`  x positions in meters with shape `[length(t) nParticles]`
+ `y`  y positions in meters with shape `[length(t) nParticles]`

## Discussion

  This method removes any initial positions that lie outside the model
  domain or inside model obstacles before integrating. The trajectories are
  sampled on the uniform output grid `0:dt:T`.
