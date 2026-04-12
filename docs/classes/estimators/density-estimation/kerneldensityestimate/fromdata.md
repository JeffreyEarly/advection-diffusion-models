---
layout: default
title: fromData
parent: KernelDensityEstimate
grand_parent: Classes
nav_order: 7
mathjax: true
---

#  fromData

Fit a kernel density estimate from sampled data.


---

## Declaration
```matlab
 self = KernelDensityEstimate.fromData(data,minimum=...,maximum=...,bandwidthGridSize=...)
```
## Parameters
+ `data`  sample locations stored as a vector for one dimension or an `N-by-2` matrix for two dimensions
+ `minimum`  optional lower endpoint of the default evaluation box; the default is `min(data) - range(data)/2` componentwise
+ `maximum`  optional upper endpoint of the default evaluation box; the default is `max(data) + range(data)/2` componentwise
+ `bandwidthGridSize`  optional DCT grid size used internally by the Botev solver; values are rounded up internally to the next power of two

## Returns
+ `self`  fitted `KernelDensityEstimate` instance

## Discussion

  `fromData(...)` chooses a Gaussian product-kernel bandwidth
  using the diffusion estimator of Botev, Grotowski, and
  Kroese, then stores the fitted data cloud together with the
  resulting bandwidth and default evaluation box.
