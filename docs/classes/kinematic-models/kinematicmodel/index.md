---
layout: default
title: KinematicModel
has_children: false
has_toc: false
mathjax: true
parent: Kinematic models
grand_parent: Class documentation
nav_order: 1
---

#  KinematicModel

Base class for two-dimensional deterministic velocity fields.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef KinematicModel < handle</code></pre></div></div>

## Overview

A kinematic model defines particle motion through

$$ \frac{dx}{dt} = u(t,x,y), \qquad \frac{dy}{dt} = v(t,x,y). $$

Subclasses specify the velocity field itself and, optionally, finite
domain limits, periodic directions, polygonal obstacles, and plotting
limits. The base class provides common geometry filtering and
visualization helpers.

```matlab
model = KinematicModel();
integrator = AdvectionDiffusionIntegrator(model, 20);
x0 = [0; 5e3];
y0 = [0; 5e3];
[~, x, y] = integrator.particleTrajectories(x0, y0, 3600, 300);
figure
model.plotTrajectories(x, y)
```




## Topics
+ Configure model domains
  + [`isSpherical`](/advection-diffusion-models/classes/kinematic-models/kinematicmodel/isspherical.html) Indicates that the x-coordinate represents longitude on a sphere.
  + [`isXPeriodic`](/advection-diffusion-models/classes/kinematic-models/kinematicmodel/isxperiodic.html) Indicates whether the x-direction wraps periodically across `xlim`.
  + [`isYPeriodic`](/advection-diffusion-models/classes/kinematic-models/kinematicmodel/isyperiodic.html) Indicates whether the y-direction wraps periodically across `ylim`.
  + [`obstacles`](/advection-diffusion-models/classes/kinematic-models/kinematicmodel/obstacles.html) Polygonal obstacle geometry as a `polyshape` array.
  + [`xlim`](/advection-diffusion-models/classes/kinematic-models/kinematicmodel/xlim.html) Finite or infinite x-domain limits in meters.
  + [`ylim`](/advection-diffusion-models/classes/kinematic-models/kinematicmodel/ylim.html) Finite or infinite y-domain limits in meters.
+ Inspect model settings
  + [`name`](/advection-diffusion-models/classes/kinematic-models/kinematicmodel/name.html) Descriptive model name used in figures and generated docs.
  + [`visualScale`](/advection-diffusion-models/classes/kinematic-models/kinematicmodel/visualscale.html) Plotting scale factor.
  + [`xVisualLimits`](/advection-diffusion-models/classes/kinematic-models/kinematicmodel/xvisuallimits.html) Finite x-limits in meters used for plotting.
  + [`yVisualLimits`](/advection-diffusion-models/classes/kinematic-models/kinematicmodel/yvisuallimits.html) Finite y-limits in meters used for plotting.
+ Evaluate velocity fields
  + [`u`](/advection-diffusion-models/classes/kinematic-models/kinematicmodel/u.html) Evaluate the x-velocity component.
  + [`v`](/advection-diffusion-models/classes/kinematic-models/kinematicmodel/v.html) Evaluate the y-velocity component.
+ Plot model diagnostics
  + [`plotBounds`](/advection-diffusion-models/classes/kinematic-models/kinematicmodel/plotbounds.html) Plot domain boundaries on the current axes.
  + [`plotTrajectories`](/advection-diffusion-models/classes/kinematic-models/kinematicmodel/plottrajectories.html) Plot particle trajectories on the current axes.
  + [`plotVelocityField`](/advection-diffusion-models/classes/kinematic-models/kinematicmodel/plotvelocityfield.html) Plot the model velocity field with quivers.
+ Work with particle domains
  + [`boundsFromLimits`](/advection-diffusion-models/classes/kinematic-models/kinematicmodel/boundsfromlimits.html) Convert rectangular limits to polygon vertices.
  + [`removeOutOfBoundsParticles`](/advection-diffusion-models/classes/kinematic-models/kinematicmodel/removeoutofboundsparticles.html) Remove initial particles outside the valid domain.


---