---
layout: default
title: GriddedStreamfunctionBootstrap
parent: GriddedStreamfunctionBootstrap
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  GriddedStreamfunctionBootstrap

Create a bootstrap from canonical restart-state properties.


---

## Declaration
```matlab
 self = GriddedStreamfunctionBootstrap(options)
```
## Parameters
+ `options.fullFit`  full-data `GriddedStreamfunction` fit
+ `options.observedTrajectories`  original drifter trajectories used for resampling
+ `options.bootstrapIndices`  whole-drifter resampling indices for each replicate
+ `options.queryTimes`  stored summary times
+ `options.scoreTimes`  stored consensus-score times
+ `options.fullSummaryUCenter`  persisted full-fit `uCenter` summary values
+ `options.fullSummaryVCenter`  persisted full-fit `vCenter` summary values
+ `options.fullSummarySigmaN`  persisted full-fit `sigma_n` summary values
+ `options.fullSummarySigmaS`  persisted full-fit `sigma_s` summary values
+ `options.fullSummaryZeta`  persisted full-fit `zeta` summary values
+ `options.summaryUCenter`  persisted bootstrap `uCenter` summary matrix
+ `options.summaryVCenter`  persisted bootstrap `vCenter` summary matrix
+ `options.summarySigmaN`  persisted bootstrap `sigma_n` summary matrix
+ `options.summarySigmaS`  persisted bootstrap `sigma_s` summary matrix
+ `options.summaryZeta`  persisted bootstrap `zeta` summary matrix
+ `options.fullFitKappaValue`  persisted full-fit scalar diffusivity diagnostic
+ `options.fullFitCoherenceValue`  persisted full-fit scalar mean coherence
+ `options.fullFitCoherenceFrequency`  persisted full-fit coherence-spectrum frequencies
+ `options.fullFitCoherenceValues`  persisted full-fit coherence-spectrum values
+ `options.bestFitKappaValue`  persisted best-fit scalar diffusivity diagnostic
+ `options.bestFitCoherenceValue`  persisted best-fit scalar mean coherence
+ `options.bestFitCoherenceFrequency`  persisted best-fit coherence-spectrum frequencies
+ `options.bestFitCoherenceValues`  persisted best-fit coherence-spectrum values
+ `options.scoreUv`  persisted bootstrap velocity consensus scores
+ `options.scoreStrain`  persisted bootstrap strain consensus scores
+ `options.scoreZeta`  persisted bootstrap vorticity consensus scores
+ `options.scoreJoint`  persisted bootstrap joint consensus scores
+ `options.bootstrapFitMetadata`  persisted exact reconstruction metadata for each bootstrap replicate
+ `options.nBootstraps`  number of bootstrap replicates
+ `options.randomSeed`  random seed used for resampling

## Returns
+ `self`  canonical `GriddedStreamfunctionBootstrap` instance

## Discussion

  Use this low-level constructor when the full fit, stored
  summaries, scores, and reconstruction metadata already exist,
  for example after reading a persisted restart file. For
  ordinary bootstrap fitting from observed trajectories, use
  `GriddedStreamfunctionBootstrap.fromTrajectories(...)`.


