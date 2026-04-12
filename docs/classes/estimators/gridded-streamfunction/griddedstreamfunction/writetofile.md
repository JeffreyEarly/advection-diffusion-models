---
layout: default
title: writeToFile
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 26
mathjax: true
---

#  writeToFile

Write this fit to a NetCDF restart file.


---

## Declaration
```matlab
ncfile = writeToFile(self,path,properties=...,shouldOverwriteExisting=...,shouldAddRequiredProperties=...,attributes=...)
```
## Parameters
+ `path` path to the NetCDF restart file to write
+ `properties` optional property names to include in addition to the required restart state
+ `shouldOverwriteExisting` (optional) boolean indicating whether to overwrite an existing file at `path`. Default `false`.
+ `shouldAddRequiredProperties` (optional) boolean indicating whether to include the required restart properties automatically. Default `true`.
+ `attributes` (optional) dictionary of additional NetCDF attributes to attach to the root group

## Returns
+ `ncfile` NetCDF file handle for the written restart file

## Discussion
  Writes the canonical solved state for this fitted estimator to a NetCDF restart file so `fromFile` can reconstruct the fit later without rerunning the trajectory fit.

```matlab
fit = GriddedStreamfunction.fromTrajectories(trajectories);
ncfile = fit.writeToFile("gridded-fit.nc", shouldOverwriteExisting=true);
fitRestart = GriddedStreamfunction.fromFile("gridded-fit.nc");
```

  Pass additional property names when you want to persist optional state beyond the required restart payload.
