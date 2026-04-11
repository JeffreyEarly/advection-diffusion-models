---
layout: default
title: fromTrajectories
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 16
mathjax: true
---

#  fromTrajectories

Fit the estimator from drifter trajectory splines.


---

## Declaration
```matlab
 self = fromTrajectories(trajectories,psiKnotPoints=...,psiS=...,fastKnotPoints=...,fastS=...,mesoscaleConstraint=...)
```
## Parameters
+ `trajectories`  nonempty vector of `TrajectorySpline` drifters
+ `psiKnotPoints`  optional cell array `{qKnot, rKnot, tKnot}` for the mesoscale basis
+ `psiS`  optional mesoscale spline degree vector `[Sq Sr St]`, default `[2 2 0]`
+ `fastKnotPoints`  optional fast temporal knot vector for COM and background
+ `fastS`  optional fast temporal spline degree, default `3`
+ `mesoscaleConstraint`  optional hard mesoscale constraint `"none"`, `"zeroVorticity"`, or `"zeroStrain"`

## Returns
+ `self`  fitted `GriddedStreamfunction` estimator

## Discussion

  Use this factory with one `TrajectorySpline` per drifter.
  Positions are sampled with `trajectory.x(trajectory.t)` and
  `trajectory.y(trajectory.t)`, and observed velocities are the
  first derivatives of the same trajectory splines.

  The fast temporal basis is used for both the center-of-mass
  trajectory and the recovered common background path. The
  tensor basis defined by `psiS` and `psiKnotPoints` is used
  only for the mesoscale streamfunction in the centered frame,
  with only the additive streamfunction gauge removed. Set
  `mesoscaleConstraint` to impose a hard zero-vorticity or
  zero-strain mesoscale fit.


