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

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef GriddedStreamfunctionBootstrap < CAAnnotatedClass</code></pre></div></div>

## Overview

`GriddedStreamfunctionBootstrap` fits one deterministic
`GriddedStreamfunction` on the full drifter ensemble, then constructs
a whole-drifter bootstrap ensemble by resampling the input
`TrajectorySpline` objects with replacement. Each replicate is
summarized by COM-local mesoscale diagnostics
$$u_c(t), v_c(t), \sigma_n(t), \sigma_s(t), \zeta(t)$$ evaluated at
the fitted center-of-mass trajectory. The class also reports the
scalar diffusivity diagnostic
$$\kappa = \langle (x_{\mathrm{sm}}^2 + y_{\mathrm{sm}}^2)/(4\Delta t)\rangle,$$
together with a scalar mean coherence over the lower-frequency half
of the mean coherence spectrum and the full mean coherence spectrum
between the fitted centered-frame mesoscale and submesoscale
velocities when the jLab coherence tools are available.

The bootstrap ensemble is ranked by a consensus score formed from
time-local `KernelDensityEstimate` fits to the bootstrap cloud. The
score uses a 2-D density fit for `(uCenter, vCenter)`, a 2-D
density fit for `(sigma_n, sigma_s)`, and a 1-D density fit for
`zeta`, omitting constrained blocks when `mesoscaleConstraint`
forces them to vanish.

The structural mesoscale degrees of freedom are determined only by
the resolved mesoscale basis and the selected hard constraint, so
the full fit and every reconstructed bootstrap replicate share the
same scalar `mesoscaleDegreesOfFreedom`.

```matlab
bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=100);
dof = bootstrap.mesoscaleDegreesOfFreedom;
bestFit = bootstrap.bestFit();
quantiles = bootstrap.summaryQuantiles([0.16 0.5 0.84]);
```




## Topics
+ Create a bootstrap ensemble
  + [`fromTrajectories`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fromtrajectories.html) Create a whole-drifter bootstrap ensemble for `GriddedStreamfunction`.
+ Read from file
  + [`fromFile`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fromfile.html) Read a bootstrap ensemble from a NetCDF restart file.
+ Write to file
  + [`writeToFile`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/writetofile.html) Write this instance to a NetCDF restart file.
+ Inspect bootstrap properties
  + [`bootstrapIndices`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bootstrapindices.html) Whole-drifter resampling indices for each bootstrap replicate.
  + [`mesoscaleDegreesOfFreedom`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/mesoscaledegreesoffreedom.html) Structural mesoscale degrees of freedom shared by the ensemble.
  + [`nBootstraps`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/nbootstraps.html) Number of bootstrap replicates in the ensemble.
  + [`observedTrajectories`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/observedtrajectories.html) Original drifter trajectories used to seed the bootstrap ensemble.
  + [`queryTimes`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/querytimes.html) Times used to store COM-local bootstrap summaries.
  + [`randomSeed`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/randomseed.html) Random-number seed used for whole-drifter resampling.
  + [`scoreTimes`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/scoretimes.html) Times used to compute the consensus score.
+ Inspect full-fit diagnostics
  + [`fullFit`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fullfit.html) Full-data gridded-streamfunction fit used as the reference solution.
  + [`fullFitCoherence`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fullfitcoherence.html) Full-fit scalar mean coherence between mesoscale and submesoscale velocities.
  + [`fullFitCoherenceSpectrum`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fullfitcoherencespectrum.html) Full-fit mean coherence spectrum on the common overlap grid.
  + [`fullFitKappa`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fullfitkappa.html) Full-fit scalar submesoscale diffusivity diagnostic.
  + [`fullSummary`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fullsummary.html) Full-data COM-local mesoscale summary evaluated on `queryTimes`.
+ Inspect best-fit diagnostics
  + [`bestBootstrapIndex`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bestbootstrapindex.html) Return the bootstrap index with the highest joint consensus score.
  + [`bestFit`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bestfit.html) Reconstruct the top-ranked bootstrap replicate.
  + [`bestFitCoherence`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bestfitcoherence.html) Best-bootstrap scalar mean coherence.
  + [`bestFitCoherenceSpectrum`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bestfitcoherencespectrum.html) Best-bootstrap mean coherence spectrum on the common overlap grid.
  + [`bestFitKappa`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bestfitkappa.html) Best-bootstrap scalar submesoscale diffusivity diagnostic.
+ Reconstruct bootstrap fits
  + [`bootstrapMetadata`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bootstrapmetadata.html) Exact reconstruction metadata for each bootstrap replicate.
  + [`fitForBootstrap`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fitforbootstrap.html) Reconstruct one bootstrap replicate exactly from saved metadata.
+ Summarize bootstrap uncertainty
  + [`bootstrapCoherence`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bootstrapcoherence.html) Lazy scalar coherence diagnostics for each bootstrap replicate.
  + [`bootstrapKappa`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bootstrapkappa.html) Lazy scalar diffusivity diagnostics for each bootstrap replicate.
  + [`scores`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/scores.html) Consensus-score components and joint score for each bootstrap fit.
  + [`summary`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/summary.html) Bootstrap COM-local mesoscale summaries evaluated on `queryTimes`.
  + [`summaryQuantiles`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/summaryquantiles.html) Return bootstrap quantiles for the stored summaries.


---