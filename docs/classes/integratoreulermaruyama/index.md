---
layout: default
title: IntegratorEulerMaruyama
has_children: false
has_toc: false
mathjax: true
parent: Class documentation
nav_order: 2
---

#  IntegratorEulerMaruyama

Integrate the Itô SDE $$dy = f(t,y)\,dt + g(t,y)\,dW_t$$.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>self = IntegratorEulerMaruyama(f,g,y0,dt=...)</code></pre></div></div>

## Overview

`IntegratorEulerMaruyama` advances the stochastic differential
equation

$$
dy = f(t,y)\,dt + g(t,y)\,dW_t
$$

with the Euler-Maruyama update

$$
y_{n+1} = y_n + f(t_n,y_n)\,dt + g(t_n,y_n)\sqrt{dt}\,\eta_n,
$$

where `\eta_n` is standard normal with the same shape as `y_n`.
`advanceToTime` preserves the existing first-order linear
interpolation between accepted stochastic steps.

    - Parameter f: deterministic drift function $$f(t,y)$$ with the same output shape as `y0`
- Parameter g: stochastic amplitude function $$g(t,y)$$ with the same output shape as `y0`
- Parameter y0: initial condition $$y_0$$ stored as `nParticles x nDims`
- Parameter dt: positive scalar timestep $$dt$$


## Topics
+ Integrators
  + [`IntegratorEulerMaruyama`](/advection-diffusion-models/classes/integratoreulermaruyama/integratoreulermaruyama.html) Create an Euler-Maruyama integrator for $$dy = f\,dt + g\,dW_t$$.
  + [`advanceOneStep`](/advection-diffusion-models/classes/integratoreulermaruyama/advanceonestep.html) Advance the SDE by one fixed timestep `dt`.
  + [`advanceToTime`](/advection-diffusion-models/classes/integratoreulermaruyama/advancetotime.html) Advance the SDE until the requested output time `t`.


---