---
layout: default
title: advanceToTime
parent: Integrator
grand_parent: Classes
nav_order: 3
mathjax: true
---

#  advanceToTime

Advance the integrator until the accepted time reaches `t`.


---

## Declaration
```matlab
 y = advanceToTime(self,t)
```
## Parameters
+ `t`  scalar target time in the same units as `dt`

## Returns
+ `y`  accepted state array $$y$$ after advancing to the requested time

## Discussion

  `advanceToTime` repeatedly applies the fixed-step update until
  `currentTime >= t`, then returns the accepted state. For the
  deterministic RK4 integrator there is no substep
  interpolation; the returned state is the most recent
  accepted state `y`.
