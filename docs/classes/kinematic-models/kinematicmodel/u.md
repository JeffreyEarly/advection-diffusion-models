---
layout: default
title: u
parent: KinematicModel
grand_parent: Classes
nav_order: 11
mathjax: true
---

#  u

Evaluate the x-velocity component.


---

## Declaration
```matlab
 uValue = u(self,t,x,y)
```
## Parameters
+ `t`  scalar evaluation time in seconds
+ `x`  x-coordinate array in meters
+ `y`  y-coordinate array in meters

## Returns
+ `uValue`  x-velocity in $$m s^{-1}$$ with the same shape as `x` and `y`

## Discussion

  Output shape matches the input `x` and `y` arrays.


