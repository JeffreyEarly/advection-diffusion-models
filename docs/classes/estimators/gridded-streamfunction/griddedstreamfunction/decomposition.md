---
layout: default
title: decomposition
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 8
mathjax: true
---

#  decomposition

Per-drifter decomposition trajectories in fixed and centered frames.


---

## Discussion

  Use `fit.decomposition` to inspect how the fitted estimator
  reconstructs the same drifters that were used in the fit. Each
  field is a `TrajectorySpline` column vector aligned one-for-one
  with `observedTrajectories`, so
  `decomposition.fixedFrame.mesoscale(k)` is the mesoscale
  trajectory for `observedTrajectories(k)`. Call
  `decomposeTrajectories` when the same fitted COM and mesoscale
  fields should be applied to a different set of drifters.

  In the fixed frame,
  `decomposition.fixedFrame.background`,
  `decomposition.fixedFrame.mesoscale`, and
  `decomposition.fixedFrame.submesoscale` satisfy

  $$
  x_k = x_k^{\mathrm{bg}} + x_k^{\mathrm{meso}} + x_k^{\mathrm{sm}},
  \qquad
  y_k = y_k^{\mathrm{bg}} + y_k^{\mathrm{meso}} + y_k^{\mathrm{sm}}.
  $$

  The background path is the common fitted background trajectory
  re-anchored so that
  $$x_k^{\mathrm{bg}}(t_{k,0}) = y_k^{\mathrm{bg}}(t_{k,0}) = 0$$
  at the first sample time $$t_{k,0}$$ of drifter $$k$$. The
  fixed-frame mesoscale trajectory carries the observed initial
  drifter position, and the fixed-frame submesoscale residual is
  zero-anchored.

  In the centered frame,
  `decomposition.centeredFrame.mesoscale` and
  `decomposition.centeredFrame.submesoscale` satisfy

  $$
  \tilde{x}_k = \tilde{x}_k^{\mathrm{meso}} + \tilde{x}_k^{\mathrm{sm}},
  \qquad
  \tilde{y}_k = \tilde{y}_k^{\mathrm{meso}} + \tilde{y}_k^{\mathrm{sm}},
  $$

  after subtracting `centerOfMassTrajectory` from the observed
  drifter path. Use these centered-frame splines when you want to
  compare the fitted mesoscale motion against the observed
  COM-relative motion.

  ```matlab
  iDrifter = 1;
  trajectory = fit.observedTrajectories(iDrifter);
  ti = trajectory.t;

  background = fit.decomposition.fixedFrame.background(iDrifter);
  mesoscale = fit.decomposition.fixedFrame.mesoscale(iDrifter);
  submesoscale = fit.decomposition.fixedFrame.submesoscale(iDrifter);
  xRecon = background.x(ti) + mesoscale.x(ti) + submesoscale.x(ti);

  centeredMesoscale = fit.decomposition.centeredFrame.mesoscale(iDrifter);
  centeredSubmesoscale = fit.decomposition.centeredFrame.submesoscale(iDrifter);
  [~, qObs, rObs] = fit.centeredCoordinates(ti, trajectory.x(ti), trajectory.y(ti));
  qRecon = centeredMesoscale.x(ti) + centeredSubmesoscale.x(ti);
  rRecon = centeredMesoscale.y(ti) + centeredSubmesoscale.y(ti);
  ```


