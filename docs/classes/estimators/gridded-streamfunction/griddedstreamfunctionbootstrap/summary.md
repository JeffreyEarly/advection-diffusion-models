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
  replicate `j`.

  ```matlab
  summary = bootstrap.summary;
  medianUCenter = median(summary.uCenter, 2);
  oneTimeCloud = summary.sigma_n(10, :);
  ```
