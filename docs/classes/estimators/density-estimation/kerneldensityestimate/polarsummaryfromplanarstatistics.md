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
+ `angleReference`  optional angular reference in radians used to choose the returned angle branch

## Returns
+ `polarSummary`  radius-angle summary with mode values and lower/upper uncertainty bounds

## Discussion

  `polarSummaryFromPlanarStatistics(...)` converts the planar mode and the
  selected enclosed-mass contour into one-dimensional radius and angle
  summaries relative to `statistics.originPoint`.

  If `angleReference` is supplied, the returned angular quantities are
  shifted by integer multiples of $$2\pi$$ so the mode angle lies on the
  branch nearest that reference.

  ```matlab
  statistics = KernelDensityEstimate.planarStatisticsFromData(data);
  polarSummary = KernelDensityEstimate.polarSummaryFromPlanarStatistics(statistics);
  ```
