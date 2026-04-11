---
layout: default
title: GriddedStreamfunctionBootstrap
parent: GriddedStreamfunctionBootstrap
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  GriddedStreamfunctionBootstrap

Create a whole-drifter bootstrap ensemble for `GriddedStreamfunction`.


---

## Declaration
```matlab
 self = GriddedStreamfunctionBootstrap(trajectories,nBootstraps=...,randomSeed=...,queryTimes=...,scoreTimes=...,scoreStride=...,psiKnotPoints=...,psiS=...,fastKnotPoints=...,fastS=...,mesoscaleConstraint=...)
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

  Use this constructor when the uncertainty analysis should
  reflect sensitivity to which drifters were observed, rather
  than only the residual variability within one fitted
  ensemble.

  The constructor first fits the full-data
  `GriddedStreamfunction`, then resamples the drifter
  trajectories with replacement and fits one replicate per
  bootstrap draw. Each replicate stores COM-local mesoscale
  diagnostics on `queryTimes`, a scalar submesoscale
  diffusivity diagnostic, and the resolved spline knot vectors
  required to reconstruct the replicate exactly later.


