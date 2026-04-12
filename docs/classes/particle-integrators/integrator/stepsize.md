---
layout: default
title: stepSize
parent: Integrator
grand_parent: Classes
nav_order: 7
mathjax: true
---

#  stepSize

Fixed timestep `dt` used by every accepted step.


---

## Type
+ Class: `double`
+ Size: `(1,1)`

## Declaration
```matlab
 self.stepSize
```
## Returns
+ `stepSize`  positive scalar timestep $$dt$$

## Discussion

  This property stores the scalar timestep that appears in the RK4
  update formulas documented for `Integrator`.
