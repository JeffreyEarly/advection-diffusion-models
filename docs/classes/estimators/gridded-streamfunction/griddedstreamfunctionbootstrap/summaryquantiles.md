---
layout: default
title: summaryQuantiles
parent: GriddedStreamfunctionBootstrap
grand_parent: Classes
nav_order: 18
mathjax: true
---

#  summaryQuantiles

Return bootstrap quantiles for the stored summaries.


---

## Declaration
```matlab
 quantiles = summaryQuantiles(self,probabilities)
```
## Parameters
+ `probabilities`  vector of quantile probabilities on `[0,1]`

## Returns
+ `quantiles`  struct of bootstrap quantiles for the stored summary fields

## Discussion

  The returned struct contains the time-varying summary fields
  `uCenter`, `vCenter`, `sigma_n`, `sigma_s`, and `zeta` with
  size `[numel(queryTimes) numel(probabilities)]`, together
  with the scalar field `kappaEstimate` of size
  `[1 numel(probabilities)]`.


