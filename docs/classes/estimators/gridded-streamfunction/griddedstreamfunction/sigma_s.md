---
layout: default
title: sigma_s
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 19
mathjax: true
---

#  sigma_s

Evaluate the shear strain field $$\sigma_s = \psi_{\tilde{x}\tilde{x}} - \psi_{\tilde{y}\tilde{y}}$$.


---

## Declaration
```matlab
 values = sigma_s(self,t,x,y)
```
## Parameters
+ `t`  scalar time or array matching `x` and `y`
+ `x`  x-coordinate array in meters
+ `y`  y-coordinate array in meters

## Returns
+ `values`  shear strain in $$s^{-1}$$

## Discussion

  This is a derived diagnostic computed from the solved
  mesoscale spline.

  ```matlab
  tFit = fit.fitSupportTimes;
  xCom = fit.centerOfMassTrajectory.x(tFit);
  yCom = fit.centerOfMassTrajectory.y(tFit);
  plot(tFit, fit.sigma_s(tFit, xCom, yCom))
  xlabel("t (s)")
  ylabel("\sigma_s (s^{-1})")
  ```
