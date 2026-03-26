---
layout: default
title: plotVelocityField
parent: KinematicModel
grand_parent: Classes
nav_order: 9
mathjax: true
---

#  plotVelocityField

Plot the model velocity field with quivers.


---

## Declaration
```matlab
 plotVelocityField(self,options)
```
## Parameters
+ `options.t`  scalar plotting time in seconds
+ `options.quiverScale`  scale factor forwarded to `quiver`
+ `options.numPoints`  number of points in the y-direction

## Returns
+ `X`  optional x-grid in meters when two outputs are requested
+ `Y`  optional y-grid in meters when two outputs are requested

## Discussion

  The plotting grid uses `2*numPoints` points in x and `numPoints` points
  in y to preserve the typical wide aspect ratio.


