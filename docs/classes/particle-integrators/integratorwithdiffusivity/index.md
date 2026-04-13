---
layout: default
title: IntegratorWithDiffusivity
has_children: false
has_toc: false
mathjax: true
parent: Particle integrators
grand_parent: Class documentation
nav_order: 3
---

#  IntegratorWithDiffusivity

Integrate $$dy = f(t,y)\,dt + \sqrt{2\kappa}\,dW_t$$ in a box domain.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IntegratorWithDiffusivity < Integrator</code></pre></div></div>

## Overview

`IntegratorWithDiffusivity` combines the deterministic RK4 drift
update from `Integrator` with additive stochastic diffusion,

$$
dy = f(t,y)\,dt + \sqrt{2\kappa}\,dW_t,
$$

where `kappa` may be a scalar or a `1 x nDims` vector of componentwise
diffusivities $$\kappa_i$$. The bounds `ymin` and `ymax` define
independent reflecting or unbounded box conditions along each
coordinate, using the repository's existing reflected and wrapped
increment formulas.

```matlab
f = @(t, y) zeros(size(y));
y0 = [0.25 0.25; 0.75 0.75];
integrator = IntegratorWithDiffusivity(f, y0, dt=0.1, kappa=0.01, ymin=0, ymax=1);
y = integrator.advanceOneStep();
```




## Topics
+ Create the integrator
  + [`IntegratorWithDiffusivity`](/advection-diffusion-models/classes/particle-integrators/integratorwithdiffusivity/integratorwithdiffusivity.html) Create an additive-diffusion integrator on a box domain.
+ Inspect integrator settings
  + [`kappa`](/advection-diffusion-models/classes/particle-integrators/integratorwithdiffusivity/kappa.html) Componentwise diffusivity `kappa`.
  + [`ymax`](/advection-diffusion-models/classes/particle-integrators/integratorwithdiffusivity/ymax.html) Upper box bounds `ymax`.
  + [`ymin`](/advection-diffusion-models/classes/particle-integrators/integratorwithdiffusivity/ymin.html) Lower box bounds `ymin`.


---