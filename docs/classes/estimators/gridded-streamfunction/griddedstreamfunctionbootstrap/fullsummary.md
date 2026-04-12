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
  aligned with `queryTimes`. These are the center diagnostics of
  the deterministic `fullFit`, evaluated at the fitted
  center-of-mass trajectory
  `fullFit.centerOfMassTrajectory(queryTimes)`.
