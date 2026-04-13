---
layout: default
title: KernelDensityEstimate
has_children: false
has_toc: false
mathjax: true
parent: Density estimation
grand_parent: Estimators
nav_order: 1
---

#  KernelDensityEstimate

Fit and evaluate Gaussian kernel density estimates in one or two dimensions.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef KernelDensityEstimate</code></pre></div></div>

## Overview

`KernelDensityEstimate` stores sampled data together with a Gaussian
product-kernel bandwidth selected from the diffusion-based Botev
estimator. For fitted models created with `fromData(...)`, the
bandwidth is chosen from the same diffusion fixed-point equations as
the legacy `kde` and `kde2d` utilities, while density and CDF values
are evaluated exactly from the Gaussian kernel sum rather than from a
reconstructed DCT grid.

In one dimension the fitted density is

$$ \hat{f}(x) = \frac{1}{N h}\sum_{i=1}^{N}\phi\!\left(\frac{x - x_i}{h}\right), $$

and in two dimensions the fitted density is

$$ \hat{f}(x,y) = \frac{1}{N h_x h_y}\sum_{i=1}^{N}\phi\!\left(\frac{x - x_i}{h_x}\right)\phi\!\left(\frac{y - y_i}{h_y}\right), $$

where $$\phi$$ is the standard normal density.

```matlab
data = [sigma_n(:), sigma_s(:)];
model = KernelDensityEstimate.fromData(data);
[density, gridVectors] = model.densityOnGrid(gridSize=[192 320]);
contourf(gridVectors{1}, gridVectors{2}, density.')
stats = KernelDensityEstimate.planarStatisticsFromData(data);
polar = KernelDensityEstimate.polarSummaryFromPlanarStatistics(stats);
strain = KernelDensityEstimate.strainSummaryFromPlanarStatistics(stats);
```




## Topics
+ Create a density estimate
  + [`fromData`](/advection-diffusion-models/classes/estimators/density-estimation/kerneldensityestimate/fromdata.html) Fit a kernel density estimate from sampled data.
  + [`KernelDensityEstimate`](/advection-diffusion-models/classes/estimators/density-estimation/kerneldensityestimate/kerneldensityestimate.html) Create a kernel density estimate from solved state.
+ Inspect density-estimate properties
  + [`bandwidth`](/advection-diffusion-models/classes/estimators/density-estimation/kerneldensityestimate/bandwidth.html) Gaussian kernel bandwidth.
  + [`data`](/advection-diffusion-models/classes/estimators/density-estimation/kerneldensityestimate/data.html) Sample locations used to define the fitted kernel density estimate.
  + [`maximum`](/advection-diffusion-models/classes/estimators/density-estimation/kerneldensityestimate/maximum.html) Upper bound of the default evaluation box.
  + [`minimum`](/advection-diffusion-models/classes/estimators/density-estimation/kerneldensityestimate/minimum.html) Lower bound of the default evaluation box.
  + [`numDimensions`](/advection-diffusion-models/classes/estimators/density-estimation/kerneldensityestimate/numdimensions.html) Number of coordinates per sample.
+ Evaluate the density estimate
  + [`cdfOnGrid`](/advection-diffusion-models/classes/estimators/density-estimation/kerneldensityestimate/cdfongrid.html) Evaluate the fitted one-dimensional CDF on a regular query grid.
  + [`densityAt`](/advection-diffusion-models/classes/estimators/density-estimation/kerneldensityestimate/densityat.html) Evaluate the fitted density at arbitrary query points.
  + [`densityOnGrid`](/advection-diffusion-models/classes/estimators/density-estimation/kerneldensityestimate/densityongrid.html) Evaluate the fitted density on a regular query grid.
  + Planar summary
    + [`planarStatisticsFromData`](/advection-diffusion-models/classes/estimators/density-estimation/kerneldensityestimate/planarstatisticsfromdata.html) Summarize a planar kernel density estimate from sampled data.
  + Planar rendering
    + [`plotPlanarStatistics`](/advection-diffusion-models/classes/estimators/density-estimation/kerneldensityestimate/plotplanarstatistics.html) Plot a planar KDE summary from precomputed statistics.
  + Polar reduction
    + [`polarSummaryFromPlanarStatistics`](/advection-diffusion-models/classes/estimators/density-estimation/kerneldensityestimate/polarsummaryfromplanarstatistics.html) Reduce planar KDE statistics to radius-angle uncertainty summaries.
  + Strain reduction
    + [`strainSummaryFromPlanarStatistics`](/advection-diffusion-models/classes/estimators/density-estimation/kerneldensityestimate/strainsummaryfromplanarstatistics.html) Reduce planar KDE statistics to strain magnitude-angle summaries.


---