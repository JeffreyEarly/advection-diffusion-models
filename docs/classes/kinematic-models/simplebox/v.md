---
layout: default
title: v
parent: SimpleBox
grand_parent: Classes
nav_order: 5
mathjax: true
---

#  v

Evaluate the zero y-velocity field.


---

## Declaration
```matlab
 vValue = v(self,t,x,y)
```
## Parameters
+ `t`  scalar evaluation time in seconds
+ `x`  x-coordinate array in meters
+ `y`  y-coordinate array in meters

## Returns
+ `vValue`  y-velocity in $$m s^{-1}$$ with the same shape as `y`

## Discussion

  This model sets $$v(t,x,y) = 0$$ everywhere in the box.


