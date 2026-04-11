---
layout: default
title: centeredCoordinates
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 4
mathjax: true
---

#  centeredCoordinates

Convert fixed-frame coordinates to COM-frame coordinates.


---

## Declaration
```matlab
 [tEval,q,r] = centeredCoordinates(self,t,x,y)
```
## Parameters
+ `t`  scalar time, vector matching the first dimension of `x` and `y`, or array matching `x` and `y`
+ `x`  x-coordinate array in meters
+ `y`  y-coordinate array in meters

## Returns
+ `tEval`  expanded evaluation-time array matching `q` and `r`
+ `q`  centered x-coordinate array $$\tilde{x}$$
+ `r`  centered y-coordinate array $$\tilde{y}$$

## Discussion

  This method evaluates the fitted center-of-mass trajectory and returns
  the centered coordinates
  $$\tilde{x} = x - m_x(t)$$ and $$\tilde{y} = y - m_y(t)$$ with the same
  shape rules used by the mesoscale evaluators.


