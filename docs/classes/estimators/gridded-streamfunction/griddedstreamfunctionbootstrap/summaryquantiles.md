---
layout: default
title: summaryQuantiles
parent: GriddedStreamfunctionBootstrap
grand_parent: Classes
nav_order: 26
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
  with the scalar fields `kappa` and `coherence` of size
  `[1 numel(probabilities)]`. The scalar coherence values use
  the lower-frequency half of each bootstrap replicate's mean
  coherence spectrum.

  ```matlab
  quantiles = bootstrap.summaryQuantiles([0.16 0.5 0.84]);
  medianStrain = quantiles.sigma_n(:, 2);
  kappaBand = quantiles.kappa;
  ```
