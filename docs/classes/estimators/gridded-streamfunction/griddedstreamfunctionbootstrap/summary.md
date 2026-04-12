---
layout: default
title: summary
parent: GriddedStreamfunctionBootstrap
grand_parent: Classes
nav_order: 25
mathjax: true
---

#  summary

Bootstrap COM-local mesoscale summaries evaluated on `queryTimes`.


---

## Discussion

  `summary.uCenter`, `summary.vCenter`, `summary.sigma_n`,
  `summary.sigma_s`, and `summary.zeta` are arrays of size
  `[numel(queryTimes) nBootstraps]`. Row `i` corresponds to
  `queryTimes(i)`, and column `j` corresponds to bootstrap
  replicate `j`. These are derived ensemble summaries for
  uncertainty analysis; use `fitForBootstrap` or `bestFit` when
  you need a full reconstructed fit rather than the stored summary
  arrays.

  ```matlab
  summary = bootstrap.summary;
  medianUCenter = median(summary.uCenter, 2);
  plot(bootstrap.queryTimes, medianUCenter)
  ```
