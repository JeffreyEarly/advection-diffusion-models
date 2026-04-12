---
layout: default
title: currentTime
parent: Integrator
grand_parent: Classes
nav_order: 4
mathjax: true
---

#  currentTime

Current integration time `t`.


---

## Type
+ Class: `double`
+ Size: `(1,1)`

## Declaration
```matlab
 self.currentTime
```
## Returns
+ `currentTime`  scalar time $$t$$ associated with `currentY`

## Discussion

  `currentTime` is updated after every accepted step and has the
  same physical units as the `t` inputs passed to the public
  stepping methods.
