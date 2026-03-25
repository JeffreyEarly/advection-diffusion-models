---
layout: default
title: IntegratorWithObstacles
has_children: false
has_toc: false
mathjax: true
parent: Class documentation
nav_order: 4
---

#  IntegratorWithObstacles

Integrate diffusion with reflecting polygonal obstacles in 2D.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>self = IntegratorWithObstacles(f,y0,dt=...,kappa=...,ymin=...,ymax=...,obstacles=polyshape.empty(0,1),isPeriodic=false)</code></pre></div></div>

## Overview

`IntegratorWithObstacles` advances the stochastic model

$$
dy = f(t,y)\,dt + \sqrt{2\kappa}\,dW_t
$$

on a two-dimensional domain with coordinate bounds `ymin` and `ymax`,
optional periodic wrapping, and polygonal obstacles. Proposed step
segments that cross an obstacle are handled by repeated geometric
reflection against the obstacle boundary.

    - Parameter f: deterministic drift function $$f(t,y)$$
- Parameter y0: initial condition $$y_0$$ stored as `nParticles x 2`
- Parameter dt: positive scalar timestep $$dt$$
- Parameter kappa: scalar or `1 x 2` diffusivity vector $$\kappa$$; defaults to `0`
- Parameter ymin: scalar or `1 x 2` lower coordinate bounds; defaults to `-Inf`
- Parameter ymax: scalar or `1 x 2` upper coordinate bounds; defaults to `Inf`
- Parameter obstacles: `polyshape` array of non-overlapping reflecting obstacles
- Parameter isPeriodic: logical scalar or `1 x 2` vector that marks periodic coordinate directions; defaults to `false`


## Topics
+ Integrators
  + [`IntegratorWithObstacles`](/advection-diffusion-models/classes/integratorwithobstacles/integratorwithobstacles.html) Create a 2D obstacle-aware stochastic integrator.
  + State
    + [`isPeriodic`](/advection-diffusion-models/classes/integratorwithobstacles/isperiodic.html) Periodicity flags for each coordinate direction.
    + [`kappa`](/advection-diffusion-models/classes/integratorwithobstacles/kappa.html) Componentwise diffusivity `kappa`.
    + [`obstacles`](/advection-diffusion-models/classes/integratorwithobstacles/obstacles.html) Polygonal obstacles used for reflecting boundaries.
    + [`ymax`](/advection-diffusion-models/classes/integratorwithobstacles/ymax.html) Upper coordinate bounds `ymax`.
    + [`ymin`](/advection-diffusion-models/classes/integratorwithobstacles/ymin.html) Lower coordinate bounds `ymin`.


---