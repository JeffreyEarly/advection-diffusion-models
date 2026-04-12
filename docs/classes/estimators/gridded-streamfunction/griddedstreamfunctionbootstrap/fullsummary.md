---
layout: default
title: fullSummary
parent: GriddedStreamfunctionBootstrap
grand_parent: Classes
nav_order: 17
mathjax: true
---

#  fullSummary

Full-data COM-local mesoscale summary evaluated on `queryTimes`.


---

## Discussion

  `fullSummary` contains the fields `uCenter`, `vCenter`,
  `sigma_n`, `sigma_s`, and `zeta`, each stored as a column vector
  aligned with `queryTimes`. These are the primary COM-local
  diagnostic outputs of the deterministic `fullFit`, evaluated at
  the fitted center-of-mass positions
  `(fullFit.centerOfMassTrajectory.x(queryTimes),`
  `fullFit.centerOfMassTrajectory.y(queryTimes))`.

  ```matlab
  fullSummary = bootstrap.fullSummary;
  plot(bootstrap.queryTimes, fullSummary.zeta)
  xlabel("t (s)")
  ylabel("\zeta_c (s^{-1})")
  ```
