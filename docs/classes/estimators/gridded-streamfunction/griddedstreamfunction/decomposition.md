---
layout: default
title: decomposition
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 6
mathjax: true
---

#  decomposition

Per-drifter decomposition trajectories in fixed and centered frames.


---

## Discussion

  `decomposition.fixedFrame.background`,
  `decomposition.fixedFrame.mesoscale`, and
  `decomposition.fixedFrame.submesoscale` are `TrajectorySpline`
  column vectors aligned with `observedTrajectories`. In the
  centered frame, `decomposition.centeredFrame.mesoscale` and
  `decomposition.centeredFrame.submesoscale` store the corresponding
  COM-frame trajectories.


