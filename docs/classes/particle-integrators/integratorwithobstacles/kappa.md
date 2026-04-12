---
layout: default
title: kappa
parent: IntegratorWithObstacles
grand_parent: Classes
nav_order: 3
mathjax: true
---

#  kappa

Componentwise diffusivity `kappa`.


---

## Type
+ Class: `double`

## Declaration
```matlab
 self.kappa
```
## Returns
+ `kappa`  scalar-expanded or componentwise diffusivity $$\kappa$$ with units of `y.^2 / t`

## Discussion

  `kappa` stores either the scalar diffusivity expanded across both
  state dimensions or the supplied `1 x 2` diffusivity vector.
