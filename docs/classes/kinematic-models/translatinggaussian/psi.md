---
layout: default
title: psi
parent: TranslatingGaussian
grand_parent: Classes
nav_order: 6
mathjax: true
---

#  psi

Evaluate the Gaussian streamfunction.


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

  The eddy center translates at the constant speed
  `[(cx) (cy)]`.
