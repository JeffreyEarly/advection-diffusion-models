---
layout: default
title: Integrator
has_children: false
has_toc: false
mathjax: true
parent: Class documentation
nav_order: 1
---

#  Integrator

Integrate the ODE $$\frac{dy}{dt} = f(t,y)$$ with fixed-step RK4.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>self = Integrator(f,y0,dt=...)</code></pre></div></div>

## Overview

`Integrator` advances the state `y` with the classical fourth-order
Runge-Kutta update

$$
\begin{aligned}
k_1 &= f(t_n,y_n), \\
k_2 &= f\!\left(t_n + \tfrac{1}{2}dt, y_n + \tfrac{1}{2}dt\,k_1\right), \\
k_3 &= f\!\left(t_n + \tfrac{1}{2}dt, y_n + \tfrac{1}{2}dt\,k_2\right), \\
k_4 &= f(t_n + dt, y_n + dt\,k_3), \\
y_{n+1} &= y_n + \tfrac{dt}{6}\left(k_1 + 2k_2 + 2k_3 + k_4\right).
\end{aligned}
$$

The state array `y` is stored as `nParticles x nDims`, and every
evaluation of `f(t,y)` must return an array with the same shape.

    - Parameter f: function handle that evaluates the ODE right-hand side $$f(t,y)$$ and returns an array with the same shape as `y0`
- Parameter y0: initial condition $$y_0$$ stored as `nParticles x nDims`
- Parameter dt: positive scalar timestep $$dt$$ used for every accepted step


## Topics
+ Integrators
  + [`Integrator`](/advection-diffusion-models/classes/integrator/integrator.html) Create an RK4 integrator for $$\frac{dy}{dt} = f(t,y)$$.
  + [`advanceOneStep`](/advection-diffusion-models/classes/integrator/advanceonestep.html) Advance the ODE by one fixed timestep `dt`.
  + [`advanceToTime`](/advection-diffusion-models/classes/integrator/advancetotime.html) Advance the integrator until the accepted time reaches `t`.
  + [`integrateToTime`](/advection-diffusion-models/classes/integrator/integratetotime.html) Integrate the ODE to each requested output time.
  + State
    + [`currentTime`](/advection-diffusion-models/classes/integrator/currenttime.html) Current integration time `t`.
    + [`currentY`](/advection-diffusion-models/classes/integrator/currenty.html) Current state `y`.
    + [`stepSize`](/advection-diffusion-models/classes/integrator/stepsize.html) Fixed timestep `dt` used by every accepted step.
    + [`totalIterations`](/advection-diffusion-models/classes/integrator/totaliterations.html) Number of accepted timesteps taken so far.


---