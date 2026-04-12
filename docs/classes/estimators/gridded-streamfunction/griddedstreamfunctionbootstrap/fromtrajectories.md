---
layout: default
title: fromTrajectories
parent: GriddedStreamfunctionBootstrap
grand_parent: Classes
nav_order: 12
mathjax: true
---

#  fromTrajectories

Create a whole-drifter bootstrap ensemble for `GriddedStreamfunction`.


---

## Declaration
```matlab
 self = fromTrajectories(trajectories,nBootstraps=...,randomSeed=...,queryTimes=...,scoreTimes=...,scoreStride=...,psiKnotPoints=...,psiS=...,fastKnotPoints=...,fastS=...,mesoscaleConstraint=...)
```
## Parameters
+ `trajectories`  nonempty vector of `TrajectorySpline` drifters
+ `nBootstraps`  number of whole-drifter bootstrap replicates, default `100`
+ `randomSeed`  integer random seed used for resampling, default `0`
+ `queryTimes`  optional strictly increasing query times for stored summaries
+ `scoreTimes`  optional strictly increasing subset of `queryTimes` used for consensus scoring
+ `scoreStride`  optional positive stride used to subsample `queryTimes` into `scoreTimes`, default `1`
+ `psiKnotPoints`  optional cell array `{qKnot, rKnot, tKnot}` for the mesoscale basis
+ `psiS`  optional mesoscale spline degree vector `[Sq Sr St]`, default `[2 2 0]`
+ `fastKnotPoints`  optional fast temporal knot vector for COM and background
+ `fastS`  optional fast temporal spline degree, default `3`
+ `mesoscaleConstraint`  optional hard mesoscale constraint `"none"`, `"zeroVorticity"`, or `"zeroStrain"`

## Returns
+ `self`  fitted bootstrap ensemble

## Discussion

  Use this factory when the uncertainty analysis should
  reflect sensitivity to which drifters were observed, rather
  than only the residual variability within one fitted
  ensemble.

  The factory first fits the full-data
  `GriddedStreamfunction`, then resamples the drifter
  trajectories with replacement and fits one replicate per
  bootstrap draw. Each replicate stores COM-local mesoscale
  diagnostics on `queryTimes`, while scalar bootstrap
  diagnostics are computed lazily on demand. The class also
  stores the resolved spline knot vectors required to
  reconstruct each replicate exactly later.

  ```matlab
  bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories( ...
      trajectories, nBootstraps=100, randomSeed=7);
  quantiles = bootstrap.summaryQuantiles([0.16 0.5 0.84]);
  plot(bootstrap.queryTimes, quantiles.uCenter(:, 2))
  xlabel("t (s)")
  ylabel("u_c (m/s)")
  ```
