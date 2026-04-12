---
layout: default
title: uBackground
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 21
mathjax: true
---

#  uBackground

Evaluate the fitted background x-velocity.


---

## Declaration
```matlab
 values = uBackground(self,t)
```
## Parameters
+ `t`  evaluation times in seconds

## Returns
+ `values`  background x-velocity in $$m s^{-1}$$

## Discussion

  `uBackground` is a derived evaluation of the solved
  `backgroundTrajectory`, not an additional fitted state
  variable.

  ```matlab
  tFit = fit.fitSupportTimes;
  plot(tFit, fit.uBackground(tFit))
  xlabel("t (s)")
  ylabel("u^{bg} (m/s)")
  ```
