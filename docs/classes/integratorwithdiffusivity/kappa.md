---
layout: default
title: kappa
parent: IntegratorWithDiffusivity
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  kappa

Componentwise diffusivity `kappa`.


---

## Declaration
```matlab
 self.kappa
```
## Returns
+ `kappa`  scalar-expanded or componentwise diffusivity $$\kappa$$ with units of `y.^2 / t`

## Discussion

  `kappa` stores either the scalar diffusivity expanded across all
  dimensions or the supplied `1 x nDims` diffusivity vector.


