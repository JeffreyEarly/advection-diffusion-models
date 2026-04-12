---
layout: default
title: densityOnGrid
parent: KernelDensityEstimate
grand_parent: Classes
nav_order: 6
mathjax: true
---

#  densityOnGrid

Evaluate the fitted density on a regular query grid.


---

## Declaration
```matlab
 [density,gridVectors] = densityOnGrid(self,gridSize=...)
```
## Parameters
+ `gridSize`  optional number of grid points; scalar in one dimension, or scalar or `[Nx Ny]` in two dimensions

## Returns
+ `density`  exact Gaussian KDE values on the regular query grid
+ `gridVectors`  cell array of one or two column vectors defining the regular query grid

## Discussion

  `densityOnGrid(...)` uses the stored default bounds together
  with an arbitrary positive evaluation grid size. The public
  evaluation grid does not need to be a power of two.
