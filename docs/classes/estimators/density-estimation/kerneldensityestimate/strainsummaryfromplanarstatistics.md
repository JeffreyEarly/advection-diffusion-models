---
layout: default
title: strainSummaryFromPlanarStatistics
parent: KernelDensityEstimate
grand_parent: Classes
nav_order: 14
mathjax: true
---

#  strainSummaryFromPlanarStatistics

Reduce planar KDE statistics to strain magnitude-angle summaries.


---

## Declaration
```matlab
 strainSummary = strainSummaryFromPlanarStatistics(statistics,thetaReference=...)
```
## Parameters
+ `statistics`  planar KDE statistics returned by `planarStatisticsFromData(...)`
+ `thetaReference`  optional `pi`-periodic extensional-axis reference angle in radians used to choose the returned angle branch

## Returns
+ `strainSummary`  strain magnitude-angle summary with mode values, lower/upper uncertainty bounds, and a zero-compatibility flag

## Discussion

  `strainSummaryFromPlanarStatistics(...)` converts the planar KDE mode and
  selected contour in the $$(\sigma_n,\sigma_s)$$ plane into the physical
  strain variables

  $$ \sigma = \sqrt{\sigma_n^2 + \sigma_s^2}, \qquad
  \theta = \tfrac{1}{2}\operatorname{atan2}(\sigma_s,\sigma_n). $$

  If `thetaReference` is supplied, the returned angular quantities are
  shifted by integer multiples of $$\pi$$ so the mode angle lies on the
  branch nearest that reference.

  If the selected contour encloses the origin, the returned magnitude
  summary remains mode-centered but is made zero-compatible by forcing the
  lower magnitude bound to zero. In that same case the angular bounds are
  reported as undefined through `thetaBounds = [NaN NaN]`.

  ```matlab
  statistics = KernelDensityEstimate.planarStatisticsFromData(data);
  strainSummary = KernelDensityEstimate.strainSummaryFromPlanarStatistics(statistics);
  ```
