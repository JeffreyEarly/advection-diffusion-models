---
layout: default
title: totalIterations
parent: Integrator
grand_parent: Classes
nav_order: 8
mathjax: true
---

#  totalIterations

Number of accepted timesteps taken so far.


---

## Type
+ Class: `double`
+ Size: `(1,1)`

## Declaration
```matlab
 self.totalIterations
```
## Returns
+ `totalIterations`  nonnegative integer count of accepted steps

## Discussion

  `totalIterations` counts accepted fixed-size updates. It does not
  count interpolated output points requested by `integrateToTime`.
