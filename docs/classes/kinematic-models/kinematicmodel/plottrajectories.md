---
layout: default
title: plotTrajectories
parent: KinematicModel
grand_parent: Classes
nav_order: 8
mathjax: true
---

#  plotTrajectories

Plot particle trajectories on the current axes.


---

## Declaration
```matlab
 plotTrajectories(self,x,y,varargin)
```
## Parameters
+ `x`  trajectory x positions in meters with shape `[nTimes nParticles]`
+ `y`  trajectory y positions in meters with shape `[nTimes nParticles]`
+ `varargin`  forwarded line styling arguments

## Discussion

  When `isXPeriodic` is true, trajectories are redrawn in wrapped copies
  so the visible path remains continuous across the plotting window.


