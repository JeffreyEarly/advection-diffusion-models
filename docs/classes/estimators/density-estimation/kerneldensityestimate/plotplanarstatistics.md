---
layout: default
title: plotPlanarStatistics
parent: KernelDensityEstimate
grand_parent: Classes
nav_order: 12
mathjax: true
---

#  plotPlanarStatistics

Plot a planar KDE summary from precomputed statistics.


---

## Declaration
```matlab
 plotPlanarStatistics(ax,statistics,samples=...,contourMasses=...,showContours=...,showScatter=...,showReferencePoint=...,showWedge=...,scatterColor=...,scatterAlpha=...,scatterSize=...,referenceFaceColor=...,referenceEdgeColor=...,referenceSize=...,wedgeColor=...,wedgeAlpha=...,wedgeResolution=...,ringRadii=...,ringColor=...,ringLineStyle=...,ringLineWidth=...,applyEqualAxes=...,applyStoredBounds=...)
```
## Parameters
+ `ax`  target axes
+ `statistics`  planar KDE statistics returned by `planarStatisticsFromData(...)`
+ `samples`  optional `N-by-2` sample cloud to draw as a scatter layer
+ `contourMasses`  enclosed-mass targets used to draw filled contours
+ `showContours`  toggle filled-density contours
+ `showScatter`  toggle the sample cloud
+ `showReferencePoint`  toggle the stored reference point marker
+ `showWedge`  toggle the wedge defined by the stored mode radius and angle bounds
+ `scatterColor`  scatter-marker face color
+ `scatterAlpha`  scatter-marker face alpha
+ `scatterSize`  scatter-marker size
+ `referenceFaceColor`  reference-point face color
+ `referenceEdgeColor`  reference-point edge color
+ `referenceSize`  reference-point marker size
+ `wedgeColor`  wedge fill color
+ `wedgeAlpha`  wedge fill alpha
+ `wedgeResolution`  number of arc samples used to draw the wedge
+ `ringRadii`  optional radii for dotted reference rings centered on `originPoint`
+ `ringColor`  reference-ring color
+ `ringLineStyle`  reference-ring line style
+ `ringLineWidth`  reference-ring line width
+ `applyEqualAxes`  toggle `axis equal`
+ `applyStoredBounds`  toggle use of `statistics.bounds`

## Discussion

  `plotPlanarStatistics(...)` draws contours, a wedge summary, a scatter
  cloud, a reference point, and optional radial guides from the data
  stored in `statistics`. Rendering geometry such as contour levels and the
  wedge polygon is built inside this helper rather than stored in the
  statistical summary.

  ```matlab
  statistics = KernelDensityEstimate.planarStatisticsFromData(data, referencePoint=[0 0]);
  KernelDensityEstimate.plotPlanarStatistics(gca, statistics, samples=data, showWedge=true);
  ```
