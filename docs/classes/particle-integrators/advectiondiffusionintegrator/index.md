---
layout: default
title: AdvectionDiffusionIntegrator
has_children: false
has_toc: false
mathjax: true
parent: Particle integrators
grand_parent: Class documentation
nav_order: 5
---

#  AdvectionDiffusionIntegrator

Integrate particles through a kinematic flow with diffusivity.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef AdvectionDiffusionIntegrator</code></pre></div></div>

## Overview

This class advances particle positions according to

$$ d\mathbf{r}(t) = \mathbf{u}(t,\mathbf{r})\,dt + \sqrt{2\kappa}\,d\mathbf{W}_t, $$

where $$\mathbf{r}(t) = [x(t)\ y(t)]$$,
$$\mathbf{u}(t,\mathbf{r})$$ is supplied by a `KinematicModel`, and
$$\kappa$$ is the scalar tracer diffusivity in $$m^2 s^{-1}$$. Boundary
limits, periodic wrapping, and polygonal obstacles are delegated to
`IntegratorWithObstacles`.

The output trajectory arrays follow the contract

- `size(t) = [nTimes 1]`
- `size(x) = [nTimes nParticles]`
- `size(y) = [nTimes nParticles]`

```matlab
model = SimpleBox();
integrator = AdvectionDiffusionIntegrator(model, 20);
x0 = [0.25; 0.75] * model.Lx;
y0 = [0.25; 0.75] * model.Ly;
[~, x, y] = integrator.particleTrajectories(x0, y0, 6 * 3600, 300);
figure
model.plotTrajectories(x, y)
```




## Topics
+ Create particle integrators
  + [`AdvectionDiffusionIntegrator`](/advection-diffusion-models/classes/particle-integrators/advectiondiffusionintegrator/advectiondiffusionintegrator.html) Create an advection-diffusion integrator.
+ Inspect integrator settings
  + [`kappa`](/advection-diffusion-models/classes/particle-integrators/advectiondiffusionintegrator/kappa.html) Scalar diffusivity in $$m^2 s^{-1}$$.
  + [`kinematicModel`](/advection-diffusion-models/classes/particle-integrators/advectiondiffusionintegrator/kinematicmodel.html) Underlying deterministic velocity model.
  + [`stepSize`](/advection-diffusion-models/classes/particle-integrators/advectiondiffusionintegrator/stepsize.html) Internal integrator step size in seconds.
+ Integrate particle trajectories
  + [`particleTrajectories`](/advection-diffusion-models/classes/particle-integrators/advectiondiffusionintegrator/particletrajectories.html) Integrate particles from initial positions.
  + [`particleTrajectoriesWithStochasticVelocity`](/advection-diffusion-models/classes/particle-integrators/advectiondiffusionintegrator/particletrajectorieswithstochasticvelocity.html) Integrate particles with extra stochastic velocity terms.


---