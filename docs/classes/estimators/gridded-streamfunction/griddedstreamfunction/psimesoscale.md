---
layout: default
title: psiMesoscale
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 15
mathjax: true
---

#  psiMesoscale

Evaluate the fitted mesoscale streamfunction.


---

## Declaration
```matlab
 psiValue = psiMesoscale(self,t,x,y)
```
## Parameters
+ `t`  scalar time or array matching `x` and `y`
+ `x`  x-coordinate array in meters
+ `y`  y-coordinate array in meters

## Returns
+ `psiValue`  mesoscale streamfunction values with the same shape as `x`

## Discussion
