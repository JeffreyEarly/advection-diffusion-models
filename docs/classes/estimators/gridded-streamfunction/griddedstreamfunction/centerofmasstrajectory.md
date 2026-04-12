---
layout: default
title: centerOfMassTrajectory
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  centerOfMassTrajectory

Fitted center-of-mass trajectory.


---

## Discussion

  This is one of the estimator's primary solved outputs.
  `centerOfMassTrajectory.x(t)` evaluates $$m_x(t)$$ and
  `centerOfMassTrajectory.y(t)` evaluates $$m_y(t)$$ on the fit
  support interval.

  ```matlab
  tFit = fit.fitSupportTimes;
  plot(fit.centerOfMassTrajectory.x(tFit), fit.centerOfMassTrajectory.y(tFit))
  axis equal
  xlabel("x (m)")
  ylabel("y (m)")
  ```
