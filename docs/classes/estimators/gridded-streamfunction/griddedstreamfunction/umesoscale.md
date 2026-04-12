---
layout: default
title: uMesoscale
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 22
mathjax: true
---

#  uMesoscale

Evaluate the mesoscale x-velocity $$-\psi_{\tilde{y}}$$.


---

## Declaration
```matlab
 uValue = uMesoscale(self,t,x,y)
```
## Parameters
+ `t`  scalar time or array matching `x` and `y`
+ `x`  x-coordinate array in meters
+ `y`  y-coordinate array in meters

## Returns
+ `uValue`  mesoscale x-velocity in $$m s^{-1}$$

## Discussion

  This is a derived velocity evaluation of the solved
  mesoscale spline in centered coordinates.

  ```matlab
  trajectory = fit.observedTrajectories(1);
  ti = trajectory.t;
  uMeso = fit.uMesoscale(ti, trajectory.x(ti), trajectory.y(ti));
  plot(ti, uMeso)
  xlabel("t (s)")
  ylabel("u^{meso} (m/s)")
  ```
