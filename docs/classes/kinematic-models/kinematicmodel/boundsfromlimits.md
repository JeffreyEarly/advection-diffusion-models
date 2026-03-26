---
layout: default
title: boundsFromLimits
parent: KinematicModel
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  boundsFromLimits

Convert rectangular limits to polygon vertices.


---

## Declaration
```matlab
 bounds = boundsFromLimits(xlim,ylim)
```
## Parameters
+ `xlim`  two-element vector of x limits in meters
+ `ylim`  two-element vector of y limits in meters

## Returns
+ `bounds`  structure with fields `xv` and `yv`

## Discussion

  Use this helper when a rectangular plotting or integration
  region should be represented as a closed polygon.


