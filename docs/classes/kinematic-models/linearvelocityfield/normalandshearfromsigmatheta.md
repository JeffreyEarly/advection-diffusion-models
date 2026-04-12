---
layout: default
title: normalAndShearFromSigmaTheta
parent: LinearVelocityField
grand_parent: Classes
nav_order: 3
mathjax: true
---

#  normalAndShearFromSigmaTheta

Convert strain magnitude and orientation to strain components.


---

## Declaration
```matlab
 [sigma_n, sigma_s] = normalAndShearFromSigmaTheta(sigma,theta)
```
## Parameters
+ `sigma`  strain magnitude in $$s^{-1}$$
+ `theta`  strain orientation in radians

## Returns
+ `sigma_n`  normal strain component in $$s^{-1}$$
+ `sigma_s`  shear strain component in $$s^{-1}$$

## Discussion

  This helper evaluates $$\sigma_n = \sigma \cos(2\theta)$$ and
  $$\sigma_s = \sigma \sin(2\theta)$$.
