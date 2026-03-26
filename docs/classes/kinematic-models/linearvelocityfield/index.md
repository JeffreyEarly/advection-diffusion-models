---
layout: default
title: LinearVelocityField
has_children: false
has_toc: false
mathjax: true
parent: Kinematic models
grand_parent: Class documentation
nav_order: 5
---

#  LinearVelocityField

Two-dimensional affine velocity field with analytical solutions.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef LinearVelocityField < StreamfunctionModel</code></pre></div></div>

## Overview

This model defines the affine flow

$$ \dot{\mathbf{r}} = \mathbf{u}_0 + A\mathbf{r}, $$

with

$$ \mathbf{u}_0 = \begin{bmatrix} u_0 \\ v_0 \end{bmatrix}, \qquad
A = \frac{1}{2}
\begin{bmatrix}
\sigma_n & \sigma_s - \zeta \\
\sigma_s + \zeta & -\sigma_n
\end{bmatrix}. $$

The strain components satisfy

$$ \sigma_n = \sigma \cos(2\theta), \qquad \sigma_s = \sigma \sin(2\theta). $$

The method `momentTensorEvolution` integrates the corresponding
second-moment system

$$ \dot{M} = AM + MA^{\top} + 2\kappa I. $$

```matlab
model = LinearVelocityField(sigma=1e-5, theta=pi / 8, zeta=0);
integrator = AdvectionDiffusionIntegrator(model, 0);
x0 = [-20e3; 20e3];
y0 = [-10e3; 10e3];
[~, x, y] = integrator.particleTrajectories(x0, y0, 12 * 3600, 600);
figure
model.plotTrajectories(x, y)
```




## Topics
+ Create the model
  + [`LinearVelocityField`](/advection-diffusion-models/classes/kinematic-models/linearvelocityfield/linearvelocityfield.html) Create a linear velocity field.
+ Inspect model parameters
  + [`sigma`](/advection-diffusion-models/classes/kinematic-models/linearvelocityfield/sigma.html) Strain magnitude in $$s^{-1}$$.
  + [`sigma_n`](/advection-diffusion-models/classes/kinematic-models/linearvelocityfield/sigma_n.html) Normal strain component $$\sigma_n = \sigma \cos(2\theta)$$.
  + [`sigma_s`](/advection-diffusion-models/classes/kinematic-models/linearvelocityfield/sigma_s.html) Shear strain component $$\sigma_s = \sigma \sin(2\theta)$$.
  + [`theta`](/advection-diffusion-models/classes/kinematic-models/linearvelocityfield/theta.html) Strain orientation in radians.
  + [`u0`](/advection-diffusion-models/classes/kinematic-models/linearvelocityfield/u0.html) Uniform background x-velocity in $$m s^{-1}$$.
  + [`v0`](/advection-diffusion-models/classes/kinematic-models/linearvelocityfield/v0.html) Uniform background y-velocity in $$m s^{-1}$$.
  + [`zeta`](/advection-diffusion-models/classes/kinematic-models/linearvelocityfield/zeta.html) Relative vorticity in $$s^{-1}$$.
+ Evaluate the velocity field
  + [`u`](/advection-diffusion-models/classes/kinematic-models/linearvelocityfield/u.html) Evaluate the affine x-velocity.
  + [`v`](/advection-diffusion-models/classes/kinematic-models/linearvelocityfield/v.html) Evaluate the affine y-velocity.
+ Evaluate the streamfunction
  + [`psi`](/advection-diffusion-models/classes/kinematic-models/linearvelocityfield/psi.html) Evaluate the quadratic streamfunction.
+ Analyze particle and moment evolution
  + [`momentTensorEvolution`](/advection-diffusion-models/classes/kinematic-models/linearvelocityfield/momenttensorevolution.html) Evolve the second-moment tensor for the linear model.
  + [`particlePath`](/advection-diffusion-models/classes/kinematic-models/linearvelocityfield/particlepath.html) Evaluate the analytical particle trajectory solution.
+ Convert flow parameters
  + [`normalAndShearFromSigmaTheta`](/advection-diffusion-models/classes/kinematic-models/linearvelocityfield/normalandshearfromsigmatheta.html) Convert strain magnitude and orientation to strain components.
  + [`sigmaThetaFromNormalAndShear`](/advection-diffusion-models/classes/kinematic-models/linearvelocityfield/sigmathetafromnormalandshear.html) Convert normal and shear strain to magnitude and orientation.


---