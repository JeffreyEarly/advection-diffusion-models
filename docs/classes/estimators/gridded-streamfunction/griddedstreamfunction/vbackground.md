---
layout: default
title: vBackground
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 23
mathjax: true
---

#  vBackground

Evaluate the fitted background y-velocity.


---

## Declaration
```matlab
 values = vBackground(self,t)
```
## Parameters
+ `t`  evaluation times in seconds

## Returns
+ `values`  background y-velocity in $$m s^{-1}$$

## Discussion

  `vBackground` is a derived evaluation of the solved
  `backgroundTrajectory`, not an additional fitted state
  variable.

  ```matlab
  tFit = fit.fitSupportTimes;
  plot(tFit, fit.vBackground(tFit))
  xlabel("t (s)")
  ylabel("v^{bg} (m/s)")
  ```
