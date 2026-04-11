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
together with a scalar mean coherence and mean coherence spectrum
between the fitted centered-frame mesoscale and submesoscale
velocities when the jLab coherence tools are available.

The bootstrap ensemble is ranked by a consensus score formed from
time-local `KernelDensityEstimate` fits to the bootstrap cloud. The
score uses a 2-D density fit for `(uCenter, vCenter)`, a 2-D
density fit for `(sigma_n, sigma_s)`, and a 1-D density fit for
`zeta`, omitting constrained blocks when `mesoscaleConstraint`
forces them to vanish.

```matlab
bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories(trajectories, nBootstraps=100);
bestFit = bootstrap.bestFit();
quantiles = bootstrap.summaryQuantiles([0.16 0.5 0.84]);
```




## Topics
+ Create a bootstrap ensemble
  + [`GriddedStreamfunctionBootstrap`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/griddedstreamfunctionbootstrap.html) Create a bootstrap from canonical restart-state properties.
  + [`fromTrajectories`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fromtrajectories.html) Create a whole-drifter bootstrap ensemble for `GriddedStreamfunction`.
+ Read from file
  + [`fromFile`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fromfile.html) Read a bootstrap ensemble from a NetCDF restart file.
+ Inspect full-fit diagnostics
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
+ Other
  + [`bestFitCoherenceFrequency`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bestfitcoherencefrequency.html) Persisted best-fit coherence-spectrum frequencies.
  + [`bestFitCoherenceFrequencyIndex`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bestfitcoherencefrequencyindex.html) Index over the best-fit coherence-spectrum frequencies.
  + [`bestFitCoherenceValue`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bestfitcoherencevalue.html) Persisted best-fit scalar mean coherence.
  + [`bestFitCoherenceValues`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bestfitcoherencevalues.html) Persisted best-fit coherence-spectrum values.
  + [`bestFitKappaValue`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bestfitkappavalue.html) Persisted best-fit scalar diffusivity diagnostic.
  + [`bootstrapCoherenceCache`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bootstrapcoherencecache.html) Optional cached bootstrap scalar coherence diagnostics.
  + [`bootstrapFitMetadata`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bootstrapfitmetadata.html) Persisted exact reconstruction metadata for each bootstrap replicate.
  + [`bootstrapIndex`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bootstrapindex.html) Index over bootstrap replicates.
  + [`bootstrapIndices`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bootstrapindices.html) Whole-drifter resampling indices for each bootstrap replicate.
  + [`bootstrapKappaCache`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/bootstrapkappacache.html) Optional cached bootstrap scalar diffusivity diagnostics.
  + [`fullFit`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fullfit.html) Full-data gridded-streamfunction fit used as the reference solution.
  + [`fullFitCoherenceFrequency`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fullfitcoherencefrequency.html) Persisted full-fit coherence-spectrum frequencies.
  + [`fullFitCoherenceFrequencyIndex`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fullfitcoherencefrequencyindex.html) Index over the full-fit coherence-spectrum frequencies.
  + [`fullFitCoherenceValue`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fullfitcoherencevalue.html) Persisted full-fit scalar mean coherence.
  + [`fullFitCoherenceValues`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fullfitcoherencevalues.html) Persisted full-fit coherence-spectrum values.
  + [`fullFitKappaValue`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fullfitkappavalue.html) Persisted full-fit scalar diffusivity diagnostic.
  + [`fullSummarySigmaN`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fullsummarysigman.html) Persisted full-fit `sigma_n` summary values.
  + [`fullSummarySigmaS`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fullsummarysigmas.html) Persisted full-fit `sigma_s` summary values.
  + [`fullSummaryUCenter`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fullsummaryucenter.html) Persisted full-fit `uCenter` summary values.
  + [`fullSummaryVCenter`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fullsummaryvcenter.html) Persisted full-fit `vCenter` summary values.
  + [`fullSummaryZeta`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/fullsummaryzeta.html) Persisted full-fit `zeta` summary values.
  + [`nBootstraps`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/nbootstraps.html) Number of bootstrap replicates in the ensemble.
  + [`observedTrajectories`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/observedtrajectories.html) Original drifter trajectories used to seed the bootstrap ensemble.
  + [`queryTimes`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/querytimes.html) Times used to store COM-local bootstrap summaries.
  + [`randomSeed`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/randomseed.html) Random-number seed used for whole-drifter resampling.
  + [`scoreJoint`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/scorejoint.html) Persisted bootstrap joint consensus scores.
  + [`scoreStrain`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/scorestrain.html) Persisted bootstrap strain consensus scores.
  + [`scoreTimes`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/scoretimes.html) Times used to compute the consensus score.
  + [`scoreUv`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/scoreuv.html) Persisted bootstrap velocity consensus scores.
  + [`scoreZeta`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/scorezeta.html) Persisted bootstrap vorticity consensus scores.
  + [`summarySigmaN`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/summarysigman.html) Persisted bootstrap `sigma_n` summary matrix.
  + [`summarySigmaS`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/summarysigmas.html) Persisted bootstrap `sigma_s` summary matrix.
  + [`summaryUCenter`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/summaryucenter.html) Persisted bootstrap `uCenter` summary matrix.
  + [`summaryVCenter`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/summaryvcenter.html) Persisted bootstrap `vCenter` summary matrix.
  + [`summaryZeta`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/summaryzeta.html) Persisted bootstrap `zeta` summary matrix.
  + [`trajectoryIndex`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/trajectoryindex.html) Index over observed drifter trajectories.
  + [`writeToFile`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunctionbootstrap/writetofile.html) Write this instance to a NetCDF restart file.


---