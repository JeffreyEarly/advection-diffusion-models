---
layout: default
title: integrateToTime
parent: Integrator
grand_parent: Classes
nav_order: 6
mathjax: true
---

#  integrateToTime

Integrate the ODE to each requested output time.


---

## Declaration
```matlab
 [y, t] = integrateToTime(self,t)
```
## Parameters
+ `t`  vector of requested output times with the same units as `dt`

## Returns
+ `y`  state history $$y(t)$$ stored as `nParticles x nDims x nTimes`
+ `t`  vector of accepted output times corresponding to the returned states

## Discussion

  The returned array `y` has shape `nParticles x nDims x nTimes`
  and stores the accepted state associated with each requested
  time entry. The returned vector `t` records the integrator
  time after each request.


