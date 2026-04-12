---
layout: default
title: zeta
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 27
mathjax: true
---

#  zeta

Evaluate the relative-vorticity field $$\zeta = \psi_{\tilde{x}\tilde{x}} + \psi_{\tilde{y}\tilde{y}}$$.


---

## Declaration
```matlab
 values = zeta(self,t,x,y)
```
## Parameters
+ `t`  scalar time or array matching `x` and `y`
+ `x`  x-coordinate array in meters
+ `y`  y-coordinate array in meters

## Returns
+ `values`  relative vorticity in $$s^{-1}$$

## Discussion

  This is a derived diagnostic computed from the solved
  mesoscale spline.

  ```matlab
  tFit = fit.fitSupportTimes;
  xCom = fit.centerOfMassTrajectory.x(tFit);
  yCom = fit.centerOfMassTrajectory.y(tFit);
  plot(tFit, fit.zeta(tFit, xCom, yCom))
  xlabel("t (s)")
  ylabel("\zeta (s^{-1})")
  ```
