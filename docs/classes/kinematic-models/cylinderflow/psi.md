---
layout: default
title: psi
parent: CylinderFlow
grand_parent: Classes
nav_order: 4
mathjax: true
---

#  psi

Evaluate the potential-flow streamfunction.


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
  $$\psi(t,x,y) = U y \left(1 - \frac{R^2}{x^2 + y^2}\right)$$.


