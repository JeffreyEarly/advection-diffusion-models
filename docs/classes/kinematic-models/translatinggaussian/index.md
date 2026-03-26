---
layout: default
title: TranslatingGaussian
has_children: false
has_toc: false
mathjax: true
parent: Kinematic models
grand_parent: Class documentation
nav_order: 9
---

#  TranslatingGaussian

Translating Gaussian eddy streamfunction.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef TranslatingGaussian < StreamfunctionModel</code></pre></div></div>

## Overview

The streamfunction is

$$ \psi(t,x,y) = e^{1/2} U L \exp\left(-\frac{(x - c_x t)^2 + (y - c_y t)^2}{2L^2}\right), $$

which produces a translating Gaussian vortex with peak velocity scale
`U` and translation velocity `(cx, cy)`.




## Topics
+ Create the model
  + [`TranslatingGaussian`](/advection-diffusion-models/classes/kinematic-models/translatinggaussian/translatinggaussian.html) Create the default translating Gaussian model.
+ Inspect model parameters
  + [`L`](/advection-diffusion-models/classes/kinematic-models/translatinggaussian/l.html) Eddy length scale in meters.
  + [`U`](/advection-diffusion-models/classes/kinematic-models/translatinggaussian/u_.html) Peak velocity scale in $$m s^-1$$.
  + [`cx`](/advection-diffusion-models/classes/kinematic-models/translatinggaussian/cx.html) Translation speed in x in $$m s^-1$$.
  + [`cy`](/advection-diffusion-models/classes/kinematic-models/translatinggaussian/cy.html) Translation speed in y in $$m s^-1$$.
+ Evaluate the streamfunction
  + [`psi`](/advection-diffusion-models/classes/kinematic-models/translatinggaussian/psi.html) Evaluate the Gaussian streamfunction.
+ Evaluate the velocity field
  + [`u`](/advection-diffusion-models/classes/kinematic-models/translatinggaussian/u.html) Evaluate the x-velocity component.
  + [`v`](/advection-diffusion-models/classes/kinematic-models/translatinggaussian/v.html) Evaluate the y-velocity component.


---