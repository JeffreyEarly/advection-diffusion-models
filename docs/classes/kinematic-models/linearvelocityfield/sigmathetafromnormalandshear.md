---
layout: default
title: sigmaThetaFromNormalAndShear
parent: LinearVelocityField
grand_parent: Classes
nav_order: 7
mathjax: true
---

#  sigmaThetaFromNormalAndShear

Convert normal and shear strain to magnitude and orientation.


---

## Declaration
```matlab
 [sigma, theta] = sigmaThetaFromNormalAndShear(sigma_n,sigma_s)
```
## Parameters
+ `sigma_n`  normal strain component in $$s^-1$$
+ `sigma_s`  shear strain component in $$s^-1$$

## Returns
+ `sigma`  strain magnitude in $$s^-1$$
+ `theta`  strain orientation in radians

## Discussion

  This helper inverts the relations
  $$\sigma_n = \sigma \cos(2\theta)$$ and
  $$\sigma_s = \sigma \sin(2\theta)$$.


