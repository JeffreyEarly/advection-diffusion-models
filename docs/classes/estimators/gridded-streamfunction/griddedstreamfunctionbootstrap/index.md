---
layout: default
title: GriddedStreamfunctionBootstrap
has_children: false
has_toc: false
mathjax: true
parent: Gridded streamfunction
grand_parent: Estimators
nav_order: 2
---

#  GriddedStreamfunctionBootstrap

Bootstrap whole-drifter gridded-streamfunction fits and consensus scores.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef GriddedStreamfunctionBootstrap < handle</code></pre></div></div>

## Overview

`GriddedStreamfunctionBootstrap` fits one deterministic
`GriddedStreamfunction` on the full drifter ensemble, then constructs
a whole-drifter bootstrap ensemble by resampling the input
`TrajectorySpline` objects with replacement. Each replicate is
summarized by COM-local mesoscale diagnostics
$$u_c(t), v_c(t), \sigma_n(t), \sigma_s(t), \zeta(t)$$ evaluated at
the fitted center-of-mass trajectory, together with a scalar
diffusivity diagnostic
$$\kappa = \langle (x_{\mathrm{sm}}^2 + y_{\mathrm{sm}}^2)/(4\Delta t)\rangle.$$

The bootstrap ensemble is ranked by a consensus score formed from
time-local `KernelDensityEstimate` fits to the bootstrap cloud. The
score uses a 2-D density fit for `(uCenter, vCenter)`, a 2-D
density fit for `(sigma_n, sigma_s)`, and a 1-D density fit for
`zeta`, omitting constrained blocks when `mesoscaleConstraint`
forces them to vanish.

```matlab
bootstrap = GriddedStreamfunctionBootstrap(trajectories, nBootstraps=100);
bestFit = bootstrap.bestFit();
quantiles = bootstrap.summaryQuantiles([0.16 0.5 0.84]);
```




## Topics
+ Create a bootstrap ensemble
  + [`GriddedStreamfunctionBootstrap`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/griddedstreamfunctionbootstrap.html) Create a whole-drifter bootstrap ensemble for `GriddedStreamfunction`.
+ Inspect bootstrap properties
  + [`bootstrapIndices`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bootstrapindices.html) Whole-drifter resampling indices for each bootstrap replicate.
  + [`fullFit`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fullfit.html) Full-data gridded-streamfunction fit used as the reference solution.
  + [`fullScalarSummary`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fullscalarsummary.html) Full-data scalar bootstrap diagnostic summary.
  + [`fullSummary`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fullsummary.html) Full-data COM-local mesoscale summary evaluated on `queryTimes`.
  + [`nBootstraps`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/nbootstraps.html) Number of bootstrap replicates in the ensemble.
  + [`observedTrajectories`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/observedtrajectories.html) Original drifter trajectories used to seed the bootstrap ensemble.
  + [`queryTimes`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/querytimes.html) Times used to store COM-local bootstrap summaries.
  + [`randomSeed`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/randomseed.html) Random-number seed used for whole-drifter resampling.
  + [`scalarSummary`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/scalarsummary.html) Scalar diagnostic summaries for each bootstrap replicate.
  + [`scoreTimes`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/scoretimes.html) Times used to compute the consensus score.
  + [`scores`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/scores.html) Consensus-score components and joint score for each bootstrap fit.
  + [`summary`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/summary.html) Bootstrap COM-local mesoscale summaries evaluated on `queryTimes`.
+ Reconstruct bootstrap fits
  + [`bestBootstrapIndex`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bestbootstrapindex.html) Return the bootstrap index with the highest joint consensus score.
  + [`bestFit`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bestfit.html) Reconstruct the top-ranked bootstrap replicate.
  + [`bootstrapMetadata`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bootstrapmetadata.html) Exact reconstruction metadata for each bootstrap replicate.
  + [`fitForBootstrap`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fitforbootstrap.html) Reconstruct one bootstrap replicate exactly from saved metadata.
+ Summarize bootstrap uncertainty
  + [`summaryQuantiles`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/summaryquantiles.html) Return bootstrap quantiles for the stored summaries.


---