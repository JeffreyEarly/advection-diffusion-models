---
layout: default
title: psi
parent: StreamfunctionModel
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  psi

Evaluate the streamfunction.


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

  Output shape matches the input `x` and `y` arrays.
