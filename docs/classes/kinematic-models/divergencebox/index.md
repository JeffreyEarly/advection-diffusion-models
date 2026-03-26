---
layout: default
title: DivergenceBox
has_children: false
has_toc: false
mathjax: true
parent: Kinematic models
grand_parent: Class documentation
nav_order: 4
---

#  DivergenceBox

Box with alternating Gaussian convergence and divergence cells.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef DivergenceBox < KinematicModel</code></pre></div></div>

## Overview

The model defines a steady velocity field by summing Gaussian cells
centered at `(x_c, y_c)`:

$$ u(x,y) = \sum_{i,j} s_i e^{1/2} U \frac{x - x_c}{L_g} \exp\left(-\frac{r^2}{2L_g^2}\right), $$

$$ v(x,y) = \sum_{i,j} s_i e^{1/2} U \frac{y - y_c}{L_g} \exp\left(-\frac{r^2}{2L_g^2}\right), $$

where `s_i = (-1)^i` alternates the sign of neighboring cells.




## Topics
+ Create the model
  + [`DivergenceBox`](/advection-diffusion-models/classes/kinematic-models/divergencebox/divergencebox.html) Create the default divergence-box model.
+ Inspect model parameters
  + [`Lg`](/advection-diffusion-models/classes/kinematic-models/divergencebox/lg.html) Gaussian length scale in meters.
  + [`Lx`](/advection-diffusion-models/classes/kinematic-models/divergencebox/lx.html) Domain width in meters.
  + [`Ly`](/advection-diffusion-models/classes/kinematic-models/divergencebox/ly.html) Domain height in meters.
  + [`U`](/advection-diffusion-models/classes/kinematic-models/divergencebox/u_.html) Velocity scale in $$m s^-1$$.
  + [`m`](/advection-diffusion-models/classes/kinematic-models/divergencebox/m.html) Number of cells along y.
  + [`n`](/advection-diffusion-models/classes/kinematic-models/divergencebox/n.html) Number of cells along x.
+ Evaluate the velocity field
  + [`u`](/advection-diffusion-models/classes/kinematic-models/divergencebox/u.html) Evaluate the x-velocity component.
  + [`v`](/advection-diffusion-models/classes/kinematic-models/divergencebox/v.html) Evaluate the y-velocity component.


---