---
layout: default
title: advanceToTime
parent: IntegratorEulerMaruyama
grand_parent: Classes
nav_order: 3
mathjax: true
---

#  advanceToTime

Advance the SDE until the requested output time `t`.


---

## Declaration
```matlab
 y = advanceToTime(self,t)
```
## Parameters
+ `t`  scalar target time in the same units as `dt`

## Returns
+ `y`  state estimate $$y(t)$$ returned with first-order interpolation between accepted steps

## Discussion

  Accepted Euler-Maruyama steps are taken until
  `currentTime >= t`. When `t` lies between the two most recent
  accepted stochastic states, the returned value is the current
  first-order linear interpolation between them.


