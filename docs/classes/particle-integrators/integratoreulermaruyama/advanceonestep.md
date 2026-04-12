---
layout: default
title: advanceOneStep
parent: IntegratorEulerMaruyama
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  advanceOneStep

Advance the SDE by one fixed timestep `dt`.


---

## Declaration
```matlab
 y = advanceOneStep(self)
```
## Returns
+ `y`  accepted state array $$y_{n+1}$$ after one stochastic timestep

## Discussion

  This method stores the previous accepted state so that a later
  call to `advanceToTime` can apply the same linear
  interpolation used by the legacy implementation.
