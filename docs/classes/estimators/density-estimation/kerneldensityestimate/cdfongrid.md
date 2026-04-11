---
layout: default
title: cdfOnGrid
parent: KernelDensityEstimate
grand_parent: Classes
nav_order: 3
mathjax: true
---

#  cdfOnGrid

Evaluate the fitted one-dimensional CDF on a regular query grid.


---

## Declaration
```matlab
 [cdfValues,gridVectors] = cdfOnGrid(self,gridSize=...)
```
## Parameters
+ `gridSize`  optional number of grid points in the regular one-dimensional query grid

## Returns
+ `cdfValues`  exact Gaussian-mixture CDF values on the regular query grid
+ `gridVectors`  one-element cell array containing the query grid as a column vector

## Discussion

  The returned CDF is the exact Gaussian-mixture cumulative
  distribution associated with the stored sample locations and
  bandwidth.


