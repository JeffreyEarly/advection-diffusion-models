---
layout: default
title: CylinderFlow
has_children: false
has_toc: false
mathjax: true
parent: Kinematic models
grand_parent: Class documentation
nav_order: 3
---

#  CylinderFlow

Potential flow around a circular cylinder.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef CylinderFlow < StreamfunctionModel</code></pre></div></div>

## Overview

The model streamfunction is

$$ \psi(t,x,y) = U y \left(1 - \frac{R^2}{x^2 + y^2}\right), $$

which gives the classical incompressible potential flow past a
cylinder of radius `R`. The cylinder interior is represented as a
polygonal obstacle.

```matlab
model = CylinderFlow();
integrator = AdvectionDiffusionIntegrator(model, 0);
x0 = [-2; -2] * model.R;
y0 = [-0.75; 0.75] * model.R;
[~, x, y] = integrator.particleTrajectories(x0, y0, 12 * 3600, 600);
figure
model.plotTrajectories(x, y)
```




## Topics
+ Create the model
  + [`CylinderFlow`](/advection-diffusion-models/classes/kinematic-models/cylinderflow/cylinderflow.html) Create the default cylinder-flow model.
+ Inspect model parameters
  + [`R`](/advection-diffusion-models/classes/kinematic-models/cylinderflow/r.html) Cylinder radius in meters.
  + [`U`](/advection-diffusion-models/classes/kinematic-models/cylinderflow/u_.html) Far-field speed in $$m s^{-1}$$.
+ Evaluate the streamfunction
  + [`psi`](/advection-diffusion-models/classes/kinematic-models/cylinderflow/psi.html) Evaluate the potential-flow streamfunction.
+ Evaluate the velocity field
  + [`u`](/advection-diffusion-models/classes/kinematic-models/cylinderflow/u.html) Evaluate the x-velocity component.
  + [`v`](/advection-diffusion-models/classes/kinematic-models/cylinderflow/v.html) Evaluate the y-velocity component.


---