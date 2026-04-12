---
layout: default
title: advanceOneStep
parent: Integrator
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  advanceOneStep

Advance the ODE by one fixed timestep `dt`.


---

## Declaration
```matlab
 y = advanceOneStep(self)
```
## Returns
+ `y`  accepted state array $$y_{n+1}$$ after one timestep

## Discussion

  This method applies one accepted RK4 update and then returns
  the new state `y`.
