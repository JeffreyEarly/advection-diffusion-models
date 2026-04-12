---
layout: default
title: backgroundTrajectory
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  backgroundTrajectory

Fitted common background trajectory.


---

## Discussion

  This is one of the estimator's primary solved outputs and stores
  the single common anchored background path recovered by the fit.
  `backgroundTrajectory.x(t)` evaluates $$x^{\mathrm{bg}}(t)$$
  and `backgroundTrajectory.y(t)` evaluates
  $$y^{\mathrm{bg}}(t)$$, with
  $$x^{\mathrm{bg}}(t_0)=0$$ and $$y^{\mathrm{bg}}(t_0)=0$$ at
  the global fit start time. The recovered background velocity is
  obtained from `backgroundTrajectory.u(t)` and
  `backgroundTrajectory.v(t)`. This is the single shared
  background path for the fitted estimator; for drifter `k`,
  `decomposition.fixedFrame.background(k)` is the same path
  re-anchored so it starts at zero at that drifter's first sample
  time.

  ```matlab
  tFit = fit.fitSupportTimes;
  plot(fit.backgroundTrajectory.x(tFit), fit.backgroundTrajectory.y(tFit))
  axis equal
  xlabel("x^{bg} (m)")
  ylabel("y^{bg} (m)")
  ```
