---
layout: default
title: GriddedStreamfunction
has_children: false
has_toc: false
mathjax: true
parent: Gridded streamfunction
grand_parent: Estimators
nav_order: 1
---

#  GriddedStreamfunction

Fit a COM-frame streamfunction estimator and trajectory decomposition.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef GriddedStreamfunction < CAAnnotatedClass</code></pre></div></div>

## Overview

`GriddedStreamfunction` fits a mesoscale streamfunction
$$\psi(\tilde{x},\tilde{y},t)$$, a center-of-mass trajectory
$$m_x(t), m_y(t)$$, and an anchored background trajectory
$$x^{\mathrm{bg}}(t), y^{\mathrm{bg}}(t)$$ from asynchronous drifter
trajectory splines using the coupled estimator described in
`asynchronous-com-fit-rev6.tex`.

The estimator uses the fast temporal spline basis both to define
the COM smoothing operator and to represent the common background
velocity. It first fits the COM trajectory, then solves one
least-squares problem for the gauge-reduced mesoscale
streamfunction coefficients in the COM frame, and finally recovers
the background velocity from the COM-space residual in that same
fast basis. Optional hard `zeroVorticity` and `zeroStrain`
constraints are applied directly to the mesoscale coefficients.

The fitted decomposition follows

$$
\dot{x} = u^{\mathrm{meso}} + u^{\mathrm{bg}} + u^{\mathrm{sm}},
\qquad
\dot{y} = v^{\mathrm{meso}} + v^{\mathrm{bg}} + v^{\mathrm{sm}},
$$

with centered coordinates
$$\tilde{x} = x - m_x(t)$$ and $$\tilde{y} = y - m_y(t).$$

The fixed-frame trajectory decomposition is stored as per-drifter
spline trajectories satisfying

$$
x_k = x_k^{\mathrm{bg}} + x_k^{\mathrm{meso}} + x_k^{\mathrm{sm}},
\qquad
y_k = y_k^{\mathrm{bg}} + y_k^{\mathrm{meso}} + y_k^{\mathrm{sm}},
$$

where the mesoscale trajectory carries the observed initial drifter
position and the background and submesoscale trajectories are
zero-anchored at the first drifter sample. In the centered frame,
the mesoscale trajectory carries the initial centered position and
the centered submesoscale trajectory is zero-anchored.

```matlab
fit = GriddedStreamfunction.fromTrajectories(trajectories);
tFit = fit.fitSupportTimes;
com = fit.centerOfMassTrajectory;
background = fit.backgroundTrajectory;
decomposition = fit.decomposition;
uMeso = fit.uMesoscale(tFit, com.x(tFit), com.y(tFit));
```




## Topics
+ Fit the estimator
  + [`fromTrajectories`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/fromtrajectories.html) Fit the estimator from drifter trajectory splines.
+ Persist and restart fits
  + [`fromFile`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/fromfile.html) Read a fitted estimator from a NetCDF restart file.
  + [`writeToFile`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/writetofile.html) Write this fit to a NetCDF restart file.
+ Inspect fit setup and structure
  + Input trajectories
    + [`observedTrajectories`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/observedtrajectories.html) Observed drifter trajectory splines used for the fit.
  + Mesoscale basis
    + [`streamfunctionSpline`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/streamfunctionspline.html) Fitted COM-frame mesoscale streamfunction spline.
    + [`psiKnotPoints`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/psiknotpoints.html) Mesoscale tensor-product knot vectors `{qKnot, rKnot, tKnot}`.
    + [`psiS`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/psis.html) Mesoscale spline degrees `[Sq Sr St]`.
    + [`mesoscaleConstraint`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/mesoscaleconstraint.html) Hard constraint applied to the fitted mesoscale streamfunction.
    + [`mesoscaleDegreesOfFreedom`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/mesoscaledegreesoffreedom.html) Identifiable mesoscale degrees of freedom after gauge and constraints.
  + Fast temporal basis
    + [`fastKnotPoints`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/fastknotpoints.html) Fast temporal knot vector used for COM and background fits.
    + [`fastS`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/fasts.html) Fast temporal spline degree for COM and background fits.
    + [`fitSupportTimes`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/fitsupporttimes.html) Sorted unique observation times used as trajectory support.
    + [`representativeTimes`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/representativetimes.html) Representative pooled times from the stride-rule fast basis.
+ Inspect primary outputs
  + [`centerOfMassTrajectory`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/centerofmasstrajectory.html) Fitted center-of-mass trajectory.
  + [`backgroundTrajectory`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/backgroundtrajectory.html) Fitted common background trajectory.
  + [`decomposition`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/decomposition.html) Per-drifter decomposition trajectories in fixed and centered frames.
+ Apply fitted decomposition
  + [`decomposeTrajectories`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/decomposetrajectories.html) Apply the fitted decomposition to a supplied drifter ensemble.
+ Evaluate derived fields
  + Coordinate transform
    + [`centeredCoordinates`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/centeredcoordinates.html) Convert fixed-frame coordinates to COM-frame coordinates.
  + Mesoscale field evaluation
    + [`psiMesoscale`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/psimesoscale.html) Evaluate the fitted mesoscale streamfunction.
    + [`uMesoscale`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/umesoscale.html) Evaluate the mesoscale x-velocity $$-\psi_{\tilde{y}}$$.
    + [`vMesoscale`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/vmesoscale.html) Evaluate the mesoscale y-velocity $$\psi_{\tilde{x}}$$.
  + Background evaluation
    + [`uBackground`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/ubackground.html) Evaluate the fitted background x-velocity.
    + [`vBackground`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/vbackground.html) Evaluate the fitted background y-velocity.
  + Derived diagnostics
    + [`sigma_n`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/sigma_n.html) Evaluate the normal strain field $$\sigma_n = -2\psi_{\tilde{x}\tilde{y}}$$.
    + [`sigma_s`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/sigma_s.html) Evaluate the shear strain field $$\sigma_s = \psi_{\tilde{x}\tilde{x}} - \psi_{\tilde{y}\tilde{y}}$$.
    + [`zeta`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/zeta.html) Evaluate the relative-vorticity field $$\zeta = \psi_{\tilde{x}\tilde{x}} + \psi_{\tilde{y}\tilde{y}}$$.
  + Visualization helper
    + [`visualPrincipalStrainAngle`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/visualprincipalstrainangle.html) Evaluate a jump-free principal strain angle for visualization.


---