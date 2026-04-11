---
layout: default
title: sigma_n
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 23
mathjax: true
---

#  sigma_n

Evaluate the normal strain field $$\sigma_n = -2\psi_{\tilde{x}\tilde{y}}$$.


---

## Declaration
```matlab
 values = sigma_n(self,t,x,y)
```
## Parameters
+ `t`  scalar time or array matching `x` and `y`
+ `x`  x-coordinate array in meters
+ `y`  y-coordinate array in meters

## Returns
+ `values`  normal strain in $$s^{-1}$$

## Discussion


