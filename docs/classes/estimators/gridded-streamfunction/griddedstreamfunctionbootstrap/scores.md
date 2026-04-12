---
layout: default
title: scores
parent: GriddedStreamfunctionBootstrap
grand_parent: Classes
nav_order: 24
mathjax: true
---

#  scores

Consensus-score components and joint score for each bootstrap fit.


---

## Discussion

  `scores.uv`, `scores.strain`, `scores.zeta`, and `scores.joint`
  are row vectors of length `nBootstraps`. These are derived
  KDE-based ranking diagnostics for the bootstrap ensemble rather
  than additional fitted state.

  ```matlab
  scores = bootstrap.scores;
  [~, iBest] = max(scores.joint);
  plot(scores.joint, ".")
  hold on
  plot(iBest, scores.joint(iBest), "o")
  ```
