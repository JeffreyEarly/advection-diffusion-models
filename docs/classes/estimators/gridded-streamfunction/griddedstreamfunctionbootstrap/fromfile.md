---
layout: default
title: fromFile
parent: GriddedStreamfunctionBootstrap
grand_parent: Classes
nav_order: 21
mathjax: true
---

#  fromFile

Read a bootstrap ensemble from a NetCDF restart file.


---

## Declaration
```matlab
 self = fromFile(path)
```
## Parameters
+ `path`  path to the NetCDF restart file

## Returns
+ `self`  reconstructed `GriddedStreamfunctionBootstrap` ensemble

## Discussion

  `fromFile` reconstructs the canonical bootstrap state written
  by `writeToFile` without rerunning the whole-drifter
  resampling workflow.


