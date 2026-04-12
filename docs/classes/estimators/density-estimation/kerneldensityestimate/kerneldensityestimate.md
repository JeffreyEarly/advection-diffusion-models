---
layout: default
title: KernelDensityEstimate
parent: KernelDensityEstimate
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  KernelDensityEstimate

Create a kernel density estimate from solved state.


---

## Declaration
```matlab
 self = KernelDensityEstimate(data=...,bandwidth=...,minimum=...,maximum=...)
```
## Parameters
+ `data`  sample locations stored as an `N-by-1` vector or `N-by-2` matrix
+ `bandwidth`  scalar Gaussian bandwidth in one dimension or `[hx hy]` in two dimensions
+ `minimum`  lower endpoint of the default evaluation box
+ `maximum`  upper endpoint of the default evaluation box

## Returns
+ `self`  canonical `KernelDensityEstimate` instance

## Discussion

  Use this cheap canonical constructor when the sample data,
  bandwidth, and default evaluation box are already known.
  Ordinary fitting workflows should prefer `fromData(...)`.
