---
layout: default
title: removeOutOfBoundsParticles
parent: KinematicModel
grand_parent: Classes
nav_order: 10
mathjax: true
---

#  removeOutOfBoundsParticles

Remove initial particles outside the valid domain.


---

## Declaration
```matlab
 [x0, y0] = removeOutOfBoundsParticles(self,x0,y0)
```
## Parameters
+ `x0`  column vector of initial x positions in meters
+ `y0`  column vector of initial y positions in meters

## Returns
+ `x0`  filtered column vector of x positions in meters
+ `y0`  filtered column vector of y positions in meters

## Discussion

  Particles are removed when they lie outside `xlim` and `ylim` or inside
  any polygonal obstacle. Periodic directions are wrapped before obstacle
  intersection tests.


