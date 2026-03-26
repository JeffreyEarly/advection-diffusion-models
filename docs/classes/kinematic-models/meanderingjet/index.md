---
layout: default
title: MeanderingJet
has_children: false
has_toc: false
mathjax: true
parent: Kinematic models
grand_parent: Class documentation
nav_order: 6
---

#  MeanderingJet

Kinematic streamfunction model for a meandering jet.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef MeanderingJet < StreamfunctionModel</code></pre></div></div>

## Overview

This model uses the Bower meandering-jet streamfunction

$$ \psi(t,x,y) = UL\left[1 - \tanh(\gamma)\right], $$

with

$$ \gamma = \frac{y - A\cos(\theta)}{L\sqrt{1 + (kA\sin(\theta))^2}}, \qquad \theta = k(x - c_x t), \qquad k = \frac{2\pi}{L_x}. $$

The resulting velocity field is incompressible and periodic in the
x-direction over the interval `[0 2*Lx]`.

References:
  Amy Bower, "A simple kinematic mechanism for mixing fluid parcels
  across a meandering jet."

```matlab
model = MeanderingJet();
integrator = AdvectionDiffusionIntegrator(model, 0);
x0 = [0.25; 1.25] * model.Lx;
y0 = [-0.5; 0.5] * model.L;
[~, x, y] = integrator.particleTrajectories(x0, y0, 3 * 86400, 1800);
figure
model.plotTrajectories(x, y)
```




## Topics
+ Create the model
  + [`MeanderingJet`](/advection-diffusion-models/classes/kinematic-models/meanderingjet/meanderingjet.html) Create the default meandering jet model.
+ Inspect model parameters
  + [`A`](/advection-diffusion-models/classes/kinematic-models/meanderingjet/a.html) Meander amplitude in meters.
  + [`L`](/advection-diffusion-models/classes/kinematic-models/meanderingjet/l.html) Jet half-width in meters.
  + [`Lx`](/advection-diffusion-models/classes/kinematic-models/meanderingjet/lx.html) Meander wavelength in meters.
  + [`U`](/advection-diffusion-models/classes/kinematic-models/meanderingjet/u_.html) Velocity scale in $$m s^{-1}$$.
  + [`cx`](/advection-diffusion-models/classes/kinematic-models/meanderingjet/cx.html) Meander phase speed in $$m s^{-1}$$.
  + [`k`](/advection-diffusion-models/classes/kinematic-models/meanderingjet/k.html) Meander wavenumber in $$m^{-1}$$.
+ Evaluate meander coordinates
  + [`gamma`](/advection-diffusion-models/classes/kinematic-models/meanderingjet/gamma.html) Evaluate the nondimensional cross-jet coordinate.
  + [`theta`](/advection-diffusion-models/classes/kinematic-models/meanderingjet/theta.html) Evaluate the meander phase.
+ Evaluate the streamfunction
  + [`psi`](/advection-diffusion-models/classes/kinematic-models/meanderingjet/psi.html) Evaluate the meandering-jet streamfunction.
+ Evaluate the velocity field
  + [`u`](/advection-diffusion-models/classes/kinematic-models/meanderingjet/u.html) Evaluate the along-jet velocity component.
  + [`v`](/advection-diffusion-models/classes/kinematic-models/meanderingjet/v.html) Evaluate the cross-jet velocity component.


---