---
layout: default
title: densityAt
parent: KernelDensityEstimate
grand_parent: Classes
nav_order: 5
mathjax: true
---

#  densityAt

Evaluate the fitted density at arbitrary query points.


---

## Declaration
```matlab
 density = densityAt(self,queryPoints)
```
## Parameters
+ `queryPoints`  evaluation points as a numeric array in one dimension or an `Nq-by-2` array in two dimensions

## Returns
+ `density`  exact Gaussian KDE values at the requested query points

## Discussion

  In one dimension `queryPoints` may be any numeric array, and
  the returned `density` matches its shape. In two dimensions
  `queryPoints` must be an `Nq-by-2` array and the returned
  `density` is an `Nq-by-1` column vector.
