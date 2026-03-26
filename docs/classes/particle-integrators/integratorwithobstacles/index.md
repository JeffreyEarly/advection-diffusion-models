---
layout: default
title: IntegratorWithObstacles
has_children: false
has_toc: false
mathjax: true
parent: Particle integrators
grand_parent: Class documentation
nav_order: 4
---

#  IntegratorWithObstacles

Integrate diffusion with reflecting polygonal obstacles in 2D.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IntegratorWithObstacles < Integrator</code></pre></div></div>

## Overview

`IntegratorWithObstacles` advances the stochastic model

$$
dy = f(t,y)\,dt + \sqrt{2\kappa}\,dW_t
$$

on a two-dimensional domain with coordinate bounds `ymin` and `ymax`,
optional periodic wrapping, and polygonal obstacles. Proposed step
segments that cross an obstacle are handled by repeated geometric
reflection against the obstacle boundary.




## Topics
+ Create the integrator
  + [`IntegratorWithObstacles`](/advection-diffusion-models/classes/particle-integrators/integratorwithobstacles/integratorwithobstacles.html) Create a 2D obstacle-aware stochastic integrator.
+ Inspect integrator settings
  + [`isPeriodic`](/advection-diffusion-models/classes/particle-integrators/integratorwithobstacles/isperiodic.html) Periodicity flags for each coordinate direction.
  + [`kappa`](/advection-diffusion-models/classes/particle-integrators/integratorwithobstacles/kappa.html) Componentwise diffusivity `kappa`.
  + [`obstacles`](/advection-diffusion-models/classes/particle-integrators/integratorwithobstacles/obstacles.html) Polygonal obstacles used for reflecting boundaries.
  + [`ymax`](/advection-diffusion-models/classes/particle-integrators/integratorwithobstacles/ymax.html) Upper coordinate bounds `ymax`.
  + [`ymin`](/advection-diffusion-models/classes/particle-integrators/integratorwithobstacles/ymin.html) Lower coordinate bounds `ymin`.


---