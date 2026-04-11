---
layout: default
title: fitForBootstrap
parent: GriddedStreamfunctionBootstrap
grand_parent: Classes
nav_order: 20
mathjax: true
---

#  fitForBootstrap

Reconstruct one bootstrap replicate exactly from saved metadata.


---

## Declaration
```matlab
 fit = fitForBootstrap(self,iBootstrap)
```
## Parameters
+ `iBootstrap`  1-based bootstrap replicate index

## Returns
+ `fit`  reconstructed `GriddedStreamfunction` replicate

## Discussion

  Use this method to recover a full `GriddedStreamfunction`
  object for a particular bootstrap replicate without storing
  every replicate fit in memory.


