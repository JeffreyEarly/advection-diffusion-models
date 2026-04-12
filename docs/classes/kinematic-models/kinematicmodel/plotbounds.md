---
layout: default
title: plotBounds
parent: KinematicModel
grand_parent: Classes
nav_order: 7
mathjax: true
---

#  plotBounds

Plot domain boundaries on the current axes.


---

## Declaration
```matlab
 plotBounds(self,varargin)
```
## Parameters
+ `varargin`  forwarded line or rectangle styling arguments

## Discussion

  Finite rectangular domains are drawn as a closed box. Semi-infinite
  domains with finite `ylim` are drawn as vertical boundary lines.
