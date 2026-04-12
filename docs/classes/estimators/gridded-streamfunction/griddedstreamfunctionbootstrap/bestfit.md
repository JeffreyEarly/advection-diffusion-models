---
layout: default
title: bestFit
parent: GriddedStreamfunctionBootstrap
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  bestFit

Reconstruct the top-ranked bootstrap replicate.


---

## Declaration
```matlab
 fit = bestFit(self)
```
## Returns
+ `fit`  reconstructed `GriddedStreamfunction` for the best bootstrap replicate

## Discussion

  The returned fit has the same structural
  `mesoscaleDegreesOfFreedom` as `self.fullFit` and
  `self.mesoscaleDegreesOfFreedom`; only the fitted spline
  coefficients and diagnostics vary across replicates.
