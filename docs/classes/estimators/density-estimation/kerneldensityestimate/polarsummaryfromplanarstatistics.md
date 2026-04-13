---
layout: default
title: polarSummaryFromPlanarStatistics
parent: KernelDensityEstimate
grand_parent: Classes
nav_order: 13
mathjax: true
---

#  polarSummaryFromPlanarStatistics

Reduce planar KDE statistics to radius-angle uncertainty summaries.


---

## Declaration
```matlab
 polarSummary = polarSummaryFromPlanarStatistics(statistics,angleReference=...)
```
## Parameters
+ `statistics`  planar KDE statistics returned by `planarStatisticsFromData(...)`
+ `angleReference`  optional `2*pi`-periodic plane-angle reference in radians used to choose the returned angle branch

## Returns
+ `polarSummary`  radius-plane-angle summary with mode values, lower/upper uncertainty bounds, and an origin-inclusion flag

## Discussion

  `polarSummaryFromPlanarStatistics(...)` converts the planar mode and the
  selected enclosed-mass contour into one-dimensional radius and angle
  summaries relative to `statistics.originPoint`. The returned angular
  quantities are the raw plane angle
  $$\alpha = \operatorname{atan2}(y-y_0, x-x_0),$$
  not any problem-specific half-angle reduction.

  If `angleReference` is supplied, the returned angular quantities are
  shifted by integer multiples of $$2\pi$$ so the mode angle lies on the
  branch nearest that reference.

  If the selected contour encloses `statistics.originPoint`, the radius
  summary remains mode-centered but the lower radial bound is set to zero
  and the angular bounds are reported as undefined. This origin-inclusive
  condition is exposed through `containsOrigin`.

  ```matlab
  statistics = KernelDensityEstimate.planarStatisticsFromData(data);
  polarSummary = KernelDensityEstimate.polarSummaryFromPlanarStatistics(statistics);
  ```
