---
layout: default
title: SimpleBox
has_children: false
has_toc: false
mathjax: true
parent: Kinematic models
grand_parent: Class documentation
nav_order: 8
---

#  SimpleBox

Rectangular box with zero deterministic velocity.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef SimpleBox < KinematicModel</code></pre></div></div>

## Overview

This model defines

$$ u(t,x,y) = 0, \qquad v(t,x,y) = 0, $$

on a finite rectangular domain. Any motion comes from the external
integrator, typically through the diffusivity term.




## Topics
+ Create the model
  + [`SimpleBox`](/advection-diffusion-models/classes/kinematic-models/simplebox/simplebox.html) Create the default box model.
+ Inspect model parameters
  + [`Lx`](/advection-diffusion-models/classes/kinematic-models/simplebox/lx.html) Domain width in meters.
  + [`Ly`](/advection-diffusion-models/classes/kinematic-models/simplebox/ly.html) Domain height in meters.
+ Evaluate the velocity field
  + [`u`](/advection-diffusion-models/classes/kinematic-models/simplebox/u.html) Evaluate the zero x-velocity field.
  + [`v`](/advection-diffusion-models/classes/kinematic-models/simplebox/v.html) Evaluate the zero y-velocity field.


---