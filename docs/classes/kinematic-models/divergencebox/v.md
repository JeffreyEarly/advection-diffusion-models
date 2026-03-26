---
layout: default
title: v
parent: DivergenceBox
grand_parent: Classes
nav_order: 9
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
+ `vValue`  y-velocity in $$m s^{-1}$$ with the same shape as `x`

## Discussion

  Neighboring Gaussian cells alternate sign along the
  x-direction.


