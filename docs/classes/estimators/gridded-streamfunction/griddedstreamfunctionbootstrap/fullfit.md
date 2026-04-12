---
layout: default
title: fullFit
parent: GriddedStreamfunctionBootstrap
grand_parent: Classes
nav_order: 13
mathjax: true
---

#  fullFit

Full-data gridded-streamfunction fit used as the reference solution.


---

## Discussion

  `fullFit` is the deterministic `GriddedStreamfunction` fit on
  the original input drifters before any resampling is applied.
  This is one of the ensemble's primary outputs.

  ```matlab
  fullFit = bootstrap.fullFit;
  tFit = fullFit.fitSupportTimes;
  plot(fullFit.centerOfMassTrajectory.x(tFit), fullFit.centerOfMassTrajectory.y(tFit))
  axis equal
  xlabel("x (m)")
  ylabel("y (m)")
  ```
