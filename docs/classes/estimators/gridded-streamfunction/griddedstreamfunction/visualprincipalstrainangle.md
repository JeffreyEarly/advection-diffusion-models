---
layout: default
title: visualPrincipalStrainAngle
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 30
mathjax: true
---

#  visualPrincipalStrainAngle

Evaluate a jump-free principal strain angle for visualization.


---

## Declaration
```matlab
 theta = visualPrincipalStrainAngle(sigma_n,sigma_s,units=...)
```
## Parameters
+ `sigma_n`  normal strain array
+ `sigma_s`  shear strain array with the same size as `sigma_n`
+ `units`  optional output units `"degrees"` or `"radians"`, default `"degrees"`

## Returns
+ `theta`  continuity-preserving principal strain angle for visualization

## Discussion

  The principal strain angle is defined by
  $$\theta_p = \tfrac{1}{2}\operatorname{atan2}(\sigma_s,\sigma_n),$$
  with a principal representative in
  $$[-45^\circ,45^\circ)$$ or $$[-\pi/4,\pi/4).$$ For plotting, this
  helper chooses the nearest equivalent branch at each step by adding
  integer multiples of $$90^\circ$$ or $$\pi/2$$ along the first
  dimension, so the plotted angle stays continuous when it crosses the
  principal-range boundary.

  Vectors are treated as one time series and preserve their input
  orientation. Matrices and higher-dimensional arrays are processed
  independently down the first dimension.


