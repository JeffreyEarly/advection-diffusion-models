---
layout: default
title: planarStatisticsFromData
parent: KernelDensityEstimate
grand_parent: Classes
nav_order: 11
mathjax: true
---

#  planarStatisticsFromData

Summarize a planar kernel density estimate from sampled data.


---

## Declaration
```matlab
 statistics = planarStatisticsFromData(data,referencePoint=...,originPoint=...,minimum=...,maximum=...,bandwidthGridSize=...,gridSize=...,summaryMass=...,boundPaddingFactor=...,minimumHalfWidth=...)
```
## Parameters
+ `data`  sampled planar coordinates as an `N-by-2` array
+ `referencePoint`  optional planar reference point stored in the returned summary
+ `originPoint`  planar origin used for radius and angle diagnostics
+ `minimum`  optional lower evaluation bound `[xmin ymin]`
+ `maximum`  optional upper evaluation bound `[xmax ymax]`
+ `bandwidthGridSize`  optional Botev-solver DCT grid size passed to `fromData(...)`
+ `gridSize`  optional evaluation grid size passed to `densityOnGrid(...)`
+ `summaryMass`  enclosed-mass fraction used to select the summary contour
+ `boundPaddingFactor`  multiplicative padding applied to the symmetric default bounds
+ `minimumHalfWidth`  lower bound for the symmetric default half-width

## Returns
+ `statistics`  planar KDE summary with density grid, mode point, selected contour, and radial/angular bounds

## Discussion

  `planarStatisticsFromData(...)` fits a two-dimensional Gaussian KDE,
  evaluates it on a regular grid, finds the KDE mode, and extracts one
  enclosed-mass contour around that mode. The returned summary stores only
  data-derived quantities, not rendering geometry.

  Relative to `originPoint = [x_0 y_0]`, the returned polar diagnostics use

  $$ r = \sqrt{(x-x_0)^2 + (y-y_0)^2}, $$

  and

  $$ \theta = \operatorname{atan2}(y-y_0, x-x_0). $$

  When explicit bounds are omitted, the default query box is the symmetric
  box centered on `originPoint` whose half-width is
  `boundPaddingFactor * max(abs(coordinates - originPoint), [], "all")`,
  clipped below by `minimumHalfWidth`.

  ```matlab
  data = [sigma_n(:), sigma_s(:)];
  statistics = KernelDensityEstimate.planarStatisticsFromData(data, referencePoint=[0 0]);
  ```
