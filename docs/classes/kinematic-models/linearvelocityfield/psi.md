---
layout: default
title: psi
parent: LinearVelocityField
grand_parent: Classes
nav_order: 5
mathjax: true
---

#  psi

Evaluate the quadratic streamfunction.


---

## Declaration
```matlab
 psiValue = psi(self,t,x,y)
```
## Parameters
+ `t`  scalar evaluation time in seconds
+ `x`  x-coordinate array in meters
+ `y`  y-coordinate array in meters

## Returns
+ `psiValue`  streamfunction values with the same shape as `x` and `y`

## Discussion

  The implemented streamfunction is

  $$ \psi = -u_0 y + v_0 x + \frac{1}{4}(\sigma_s + \zeta)x^2 - \frac{1}{2}\sigma_n xy - \frac{1}{4}(\sigma_s - \zeta)y^2. $$


