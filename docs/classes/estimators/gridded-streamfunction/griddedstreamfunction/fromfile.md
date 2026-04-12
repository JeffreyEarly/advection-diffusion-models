---
layout: default
title: fromFile
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 9
mathjax: true
---

#  fromFile

Read a fitted estimator from a NetCDF restart file.


---

## Declaration
```matlab
 self = fromFile(path)
```
## Parameters
+ `path`  path to the NetCDF restart file

## Returns
+ `self`  reconstructed `GriddedStreamfunction` estimator

## Discussion

  `fromFile` reconstructs the canonical solved state written by
  `writeToFile` without rerunning the trajectory fit.
