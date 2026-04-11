---
layout: default
title: GriddedStreamfunction
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  GriddedStreamfunction

Create a fit from canonical solved-state properties.


---

## Declaration
```matlab
 self = GriddedStreamfunction(options)
```
## Parameters
+ `options.streamfunctionSpline`  fitted centered-frame mesoscale streamfunction spline
+ `options.observedTrajectories`  observed drifter trajectory splines aligned with the fit
+ `options.centerOfMassTrajectory`  fitted center-of-mass trajectory
+ `options.backgroundTrajectory`  fitted common background trajectory
+ `options.fixedFrameBackgroundTrajectories`  stored fixed-frame background decomposition trajectories
+ `options.fixedFrameMesoscaleTrajectories`  stored fixed-frame mesoscale decomposition trajectories
+ `options.fixedFrameSubmesoscaleTrajectories`  stored fixed-frame submesoscale decomposition trajectories
+ `options.centeredFrameMesoscaleTrajectories`  stored centered-frame mesoscale decomposition trajectories
+ `options.centeredFrameSubmesoscaleTrajectories`  stored centered-frame submesoscale decomposition trajectories
+ `options.mesoscaleConstraint`  hard mesoscale constraint `"none"`, `"zeroVorticity"`, or `"zeroStrain"`
+ `options.representativeTimes`  pooled representative times from the stride-rule fast basis
+ `options.fitSupportTimes`  sorted unique observation times used as trajectory support

## Returns
+ `self`  canonical `GriddedStreamfunction` instance

## Discussion

  Use this low-level constructor when the fitted spline state
  is already available, for example after reading a persisted
  restart file. For fitting from observed drifter trajectories,
  use `GriddedStreamfunction.fromTrajectories(...)`.

  ```matlab
  fit = GriddedStreamfunction( ...
      streamfunctionSpline=psiSpline, ...
      observedTrajectories=trajectories, ...
      centerOfMassTrajectory=centerTrajectory, ...
      backgroundTrajectory=backgroundTrajectory, ...
      fixedFrameBackgroundTrajectories=backgroundTrajectories, ...
      fixedFrameMesoscaleTrajectories=mesoscaleTrajectories, ...
      fixedFrameSubmesoscaleTrajectories=submesoscaleTrajectories, ...
      centeredFrameMesoscaleTrajectories=centeredMesoscaleTrajectories, ...
      centeredFrameSubmesoscaleTrajectories=centeredSubmesoscaleTrajectories, ...
      mesoscaleConstraint="none", ...
      representativeTimes=representativeTimes, ...
      fitSupportTimes=fitSupportTimes);
  ```


