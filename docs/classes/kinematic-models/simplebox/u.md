---
layout: default
title: u
parent: SimpleBox
grand_parent: Classes
nav_order: 4
mathjax: true
---

#  u

Evaluate the zero x-velocity field.


---

## Declaration
```matlab
 uValue = u(self,t,x,y)
```
## Parameters
+ `t`  scalar evaluation time in seconds
+ `x`  x-coordinate array in meters
+ `y`  y-coordinate array in meters

## Returns
+ `uValue`  x-velocity in $$m s^-1$$ with the same shape as `x`

## Discussion

  This model sets $$u(t,x,y) = 0$$ everywhere in the box.


