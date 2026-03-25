---
layout: default
title: IntegratorWithDiffusivity
has_children: false
has_toc: false
mathjax: true
parent: Class documentation
nav_order: 3
---

#  IntegratorWithDiffusivity

Integrate $$dy = f(t,y)\,dt + \sqrt{2\kappa}\,dW_t$$ in a box domain.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>self = IntegratorWithDiffusivity(f,y0,dt=...,kappa=...,ymin=...,ymax=...)</code></pre></div></div>

## Overview

`IntegratorWithDiffusivity` combines the deterministic RK4 drift
update from `Integrator` with additive stochastic diffusion,

$$
dy = f(t,y)\,dt + \sqrt{2\kappa}\,dW_t,
$$

where `kappa` may be a scalar or a `1 x nDims` vector of componentwise
diffusivities `\kappa_i`. The bounds `ymin` and `ymax` define
independent reflecting or unbounded box conditions along each
coordinate, using the repository's existing reflected and wrapped
increment formulas.

    - Parameter f: deterministic drift function $$f(t,y)$$
- Parameter y0: initial condition $$y_0$$ stored as `nParticles x nDims`
- Parameter dt: positive scalar timestep $$dt$$
- Parameter kappa: scalar or `1 x nDims` diffusivity vector $$\kappa$$ with units of `y.^2 / t`; defaults to `0`
- Parameter ymin: scalar or `1 x nDims` lower box bounds; defaults to `-Inf`
- Parameter ymax: scalar or `1 x nDims` upper box bounds; defaults to `Inf`


## Topics
+ Integrators
  + [`IntegratorWithDiffusivity`](/advection-diffusion-models/classes/integratorwithdiffusivity/integratorwithdiffusivity.html) Create an additive-diffusion integrator on a box domain.
  + State
    + [`kappa`](/advection-diffusion-models/classes/integratorwithdiffusivity/kappa.html) Componentwise diffusivity `kappa`.
    + [`ymax`](/advection-diffusion-models/classes/integratorwithdiffusivity/ymax.html) Upper box bounds `ymax`.
    + [`ymin`](/advection-diffusion-models/classes/integratorwithdiffusivity/ymin.html) Lower box bounds `ymin`.


---