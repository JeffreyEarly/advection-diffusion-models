---
layout: default
title: v
parent: LinearVelocityField
grand_parent: Classes
nav_order: 13
mathjax: true
---

#  v

Evaluate the affine y-velocity.


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
+ `vValue`  y-velocity in $$m s^{-1}$$ with the same shape as `x` and `y`

## Discussion

  The implemented velocity component is

  $$ v = v_0 + \frac{1}{2}\left((\sigma_s + \zeta)x - \sigma_n y\right). $$


