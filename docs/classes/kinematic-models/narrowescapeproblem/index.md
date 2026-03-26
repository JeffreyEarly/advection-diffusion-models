---
layout: default
title: NarrowEscapeProblem
has_children: false
has_toc: false
mathjax: true
parent: Kinematic models
grand_parent: Class documentation
nav_order: 7
---

#  NarrowEscapeProblem

Rectangular chamber with a narrow opening in the left wall.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef NarrowEscapeProblem < KinematicModel</code></pre></div></div>

## Overview

This model has zero deterministic velocity,

$$ u(t,x,y) = 0, \qquad v(t,x,y) = 0, $$

and represents the geometry through polygonal obstacles that leave a
gap of width `W` in the left wall. The problem is useful for testing
stochastic escape through a small opening.

```matlab
model = NarrowEscapeProblem();
integrator = AdvectionDiffusionIntegrator(model, 20);
x0 = [0.25; 0.75] * model.L;
y0 = [0.4; 0.6] * model.L;
[~, x, y] = integrator.particleTrajectories(x0, y0, 6 * 3600, 300);
figure
model.plotTrajectories(x, y)
```




## Topics
+ Create the model
  + [`NarrowEscapeProblem`](/advection-diffusion-models/classes/kinematic-models/narrowescapeproblem/narrowescapeproblem.html) Create the default narrow-escape geometry.
+ Inspect model parameters
  + [`L`](/advection-diffusion-models/classes/kinematic-models/narrowescapeproblem/l.html) Chamber length in meters.
  + [`W`](/advection-diffusion-models/classes/kinematic-models/narrowescapeproblem/w.html) Opening width in meters.
  + [`delta`](/advection-diffusion-models/classes/kinematic-models/narrowescapeproblem/delta.html) Wall thickness in meters.
+ Rebuild model geometry
  + [`regenerateObstacles`](/advection-diffusion-models/classes/kinematic-models/narrowescapeproblem/regenerateobstacles.html) Rebuild the chamber walls and escape opening.


---