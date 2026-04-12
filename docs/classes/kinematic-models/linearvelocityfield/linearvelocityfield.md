---
layout: default
title: LinearVelocityField
parent: LinearVelocityField
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  LinearVelocityField

Create a linear velocity field.


---

## Declaration
```matlab
 self = LinearVelocityField(sigma=...,theta=...,zeta=...,u0=...,v0=...)
```
## Parameters
+ `sigma`  optional strain magnitude in $$s^{-1}$$; the default is `0`
+ `theta`  optional strain orientation in radians; the default is `0`
+ `zeta`  optional relative vorticity in $$s^{-1}$$; the default is `0`
+ `u0`  optional uniform background x-velocity in $$m s^{-1}$$; the default is `0`
+ `v0`  optional uniform background y-velocity in $$m s^{-1}$$; the default is `0`

## Returns
+ `self`  `LinearVelocityField` instance

## Discussion

  Pass any subset of the model parameters by name. Omitted
  parameters default to zero.
