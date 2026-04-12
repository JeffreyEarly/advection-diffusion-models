---
layout: default
title: bootstrapCoherence
parent: GriddedStreamfunctionBootstrap
grand_parent: Classes
nav_order: 6
mathjax: true
---

#  bootstrapCoherence

Lazy scalar coherence diagnostics for each bootstrap replicate.


---

## Discussion

  `bootstrapCoherence` is a row vector of length `nBootstraps`.
  Each value uses the finite lower-frequency half of the stored
  mean coherence spectrum. When coherence is unavailable, the
  vector is filled with `NaN`.

  ```matlab
  coherence = bootstrap.bootstrapCoherence;
  plot(coherence, ".")
  xlabel("bootstrap replicate")
  ylabel("coherence")
  ```
