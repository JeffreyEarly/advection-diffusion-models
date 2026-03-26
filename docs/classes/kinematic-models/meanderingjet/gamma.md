---
layout: default
title: gamma
parent: MeanderingJet
grand_parent: Classes
nav_order: 7
mathjax: true
---

#  gamma

Evaluate the nondimensional cross-jet coordinate.


---

## Declaration
```matlab
 gammaValue = gamma(self,t,x,y)
```
## Parameters
+ `t`  scalar evaluation time in seconds
+ `x`  x-coordinate array in meters
+ `y`  y-coordinate array in meters

## Returns
+ `gammaValue`  nondimensional coordinate with the same shape as `x` and `y`

## Discussion

  `gamma` is the transformed cross-jet coordinate that appears
  inside the Bower streamfunction.


