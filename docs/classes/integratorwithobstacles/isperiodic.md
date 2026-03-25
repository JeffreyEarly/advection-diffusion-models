---
layout: default
title: isPeriodic
parent: IntegratorWithObstacles
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  isPeriodic

Periodicity flags for each coordinate direction.


---

## Declaration
```matlab
 self.isPeriodic
```
## Returns
+ `isPeriodic`  logical `1 x 2` vector of periodic-direction flags

## Discussion

  `isPeriodic` is normalized to a logical `1 x 2` row vector. A
  true entry means the corresponding coordinate is wrapped into the
  interval `[ymin(i), ymax(i))` before obstacle checks.


