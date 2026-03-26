---
layout: default
title: v
parent: KinematicModel
grand_parent: Classes
nav_order: 12
mathjax: true
---

#  v

Evaluate the y-velocity component.


---

## Declaration
```matlab
 vValue = v(self,t,x,y)
```
## Parameters
+ `t`  scalar evaluation time in seconds
+ `x`  x-coordinate array in meters
+ `y`  y-coordinate array in meters

## Returns
+ `vValue`  y-velocity in $$m s^-1$$ with the same shape as `x` and `y`

## Discussion

  Output shape matches the input `x` and `y` arrays.


