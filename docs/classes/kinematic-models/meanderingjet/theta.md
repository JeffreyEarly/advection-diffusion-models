---
layout: default
title: theta
parent: MeanderingJet
grand_parent: Classes
nav_order: 10
mathjax: true
---

#  theta

Evaluate the meander phase.


---

## Declaration
```matlab
 thetaValue = theta(self,t,x)
```
## Parameters
+ `t`  scalar evaluation time in seconds
+ `x`  x-coordinate array in meters

## Returns
+ `thetaValue`  meander phase in radians with the same shape as `x`

## Discussion

  The phase is $$\theta = k(x - c_x t)$$.
