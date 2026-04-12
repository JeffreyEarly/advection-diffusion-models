---
layout: default
title: plotStreamfunction
parent: StreamfunctionModel
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  plotStreamfunction

Plot streamfunction contours on the current axes.


---

## Declaration
```matlab
 plotStreamfunction(self,t)
```
## Parameters
+ `t`  scalar plotting time in seconds; the default is `0`

## Returns
+ `X`  optional x-grid in meters when two outputs are requested
+ `Y`  optional y-grid in meters when two outputs are requested

## Discussion

  The contour grid spans `xVisualLimits` and `yVisualLimits` and excludes
  obstacle interiors when polygonal obstacles are present.
