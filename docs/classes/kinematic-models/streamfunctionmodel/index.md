---
layout: default
title: StreamfunctionModel
has_children: false
has_toc: false
mathjax: true
parent: Kinematic models
grand_parent: Class documentation
nav_order: 2
---

#  StreamfunctionModel

Base class for incompressible two-dimensional streamfunction models.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef StreamfunctionModel < KinematicModel</code></pre></div></div>

## Overview

Streamfunction models define a scalar field $$\psi(t,x,y)$$ such that

$$ u(t,x,y) = -\frac{\partial \psi}{\partial y}, \qquad v(t,x,y) = \frac{\partial \psi}{\partial x}. $$

Subclasses implement `psi(self, t, x, y)` and inherit the common
streamfunction plotting helper.




## Topics
+ Evaluate streamfunctions
  + [`psi`](/advection-diffusion-models/classes/kinematic-models/streamfunctionmodel/psi.html) Evaluate the streamfunction.
+ Plot streamfunctions
  + [`plotStreamfunction`](/advection-diffusion-models/classes/kinematic-models/streamfunctionmodel/plotstreamfunction.html) Plot streamfunction contours on the current axes.


---