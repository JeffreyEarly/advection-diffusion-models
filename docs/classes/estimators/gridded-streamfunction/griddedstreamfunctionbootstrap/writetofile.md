---
layout: default
title: writeToFile
parent: GriddedStreamfunctionBootstrap
grand_parent: Classes
nav_order: 27
mathjax: true
---

#  writeToFile

Write this ensemble to a NetCDF restart file.


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
  Writes the canonical bootstrap state to a NetCDF restart file so `fromFile` can reconstruct the saved ensemble later without rerunning the whole-drifter bootstrap workflow.

```matlab
bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories( ...
    trajectories, nBootstraps=100);
ncfile = bootstrap.writeToFile("gridded-bootstrap.nc", shouldOverwriteExisting=true);
bootstrapRestart = GriddedStreamfunctionBootstrap.fromFile("gridded-bootstrap.nc");
```

  Pass additional property names when you want to persist optional state beyond the required restart payload.
