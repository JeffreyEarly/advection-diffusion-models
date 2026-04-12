---
layout: default
title: vMesoscale
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 24
mathjax: true
---

#  vMesoscale

Evaluate the mesoscale y-velocity $$\psi_{\tilde{x}}$$.


---

## Declaration
```matlab
 vValue = vMesoscale(self,t,x,y)
```
## Parameters
+ `t`  scalar time or array matching `x` and `y`
+ `x`  x-coordinate array in meters
+ `y`  y-coordinate array in meters

## Returns
+ `vValue`  mesoscale y-velocity in $$m s^{-1}$$

## Discussion

  This is a derived velocity evaluation of the solved
  mesoscale spline in centered coordinates.

  ```matlab
  trajectory = fit.observedTrajectories(1);
  ti = trajectory.t;
  vMeso = fit.vMesoscale(ti, trajectory.x(ti), trajectory.y(ti));
  plot(ti, vMeso)
  xlabel("t (s)")
  ylabel("v^{meso} (m/s)")
  ```
