---
layout: default
title: bootstrapKappa
parent: GriddedStreamfunctionBootstrap
grand_parent: Classes
nav_order: 8
mathjax: true
---

#  bootstrapKappa

Lazy scalar diffusivity diagnostics for each bootstrap replicate.


---

## Discussion

  `bootstrapKappa` is a row vector of length `nBootstraps`.

  ```matlab
  kappa = bootstrap.bootstrapKappa;
  plot(kappa, ".")
  xlabel("bootstrap replicate")
  ylabel("\kappa")
  ```
