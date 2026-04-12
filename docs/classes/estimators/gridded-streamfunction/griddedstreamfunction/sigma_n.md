---
layout: default
title: sigma_n
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 18
mathjax: true
---

#  sigma_n

Evaluate the normal strain field $$\sigma_n = -2\psi_{\tilde{x}\tilde{y}}$$.


---

## Declaration
```matlab
 values = sigma_n(self,t,x,y)
```
## Parameters
+ `t`  scalar time or array matching `x` and `y`
+ `x`  x-coordinate array in meters
+ `y`  y-coordinate array in meters

## Returns
+ `values`  normal strain in $$s^{-1}$$

## Discussion

  This is a derived diagnostic computed from the solved
  mesoscale spline.

  ```matlab
  tFit = fit.fitSupportTimes;
  xCom = fit.centerOfMassTrajectory.x(tFit);
  yCom = fit.centerOfMassTrajectory.y(tFit);
  plot(tFit, fit.sigma_n(tFit, xCom, yCom))
  xlabel("t (s)")
  ylabel("\sigma_n (s^{-1})")
  ```
