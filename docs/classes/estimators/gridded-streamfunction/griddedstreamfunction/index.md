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

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef GriddedStreamfunction < handle</code></pre></div></div>

## Overview

`GriddedStreamfunction` fits a mesoscale streamfunction
$$\psi(\tilde{x},\tilde{y},t)$$, a center-of-mass trajectory
$$m_x(t), m_y(t)$$, and an anchored background trajectory
$$x^{\mathrm{bg}}(t), y^{\mathrm{bg}}(t)$$ from asynchronous drifter
trajectory splines using the rev3 mesoscale-only formulation
described in `asynchronous-com-fit-rev3.tex`, with optional hard
mesoscale constraints described in `asynchronous-com-fit-rev4.tex`.

The estimator first fits the COM position in a fast temporal basis,
then solves one least-squares problem for the mesoscale
streamfunction coefficients. The background velocity is recovered
afterward as the COM residual projected back onto the same fast
basis, then integrated to obtain the common anchored background
path. When requested, the mesoscale fit imposes a hard
`zeroVorticity` or `zeroStrain` constraint on the gauge-reduced
streamfunction coefficients.

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
fit = GriddedStreamfunction(trajectories);
uMeso = fit.uMesoscale(tQuery, xQuery, yQuery);
xMeso = fit.decomposition.fixedFrame.mesoscale(1).x(tQuery);
decomposition = fit.decomposeTrajectories(otherTrajectories);
```




## Topics
+ Fit the estimator
  + [`GriddedStreamfunction`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/griddedstreamfunction.html) Fit the estimator from drifter trajectory splines.
+ Inspect fitted components
  + [`backgroundTrajectory`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/backgroundtrajectory.html) Fitted common background trajectory.
  + [`centerOfMassTrajectory`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/centerofmasstrajectory.html) Fitted center-of-mass trajectory.
  + [`fastKnotPoints`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/fastknotpoints.html) Fast temporal knot vector used for COM and background fits.
  + [`fastS`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/fasts.html) Fast temporal spline degree for COM and background fits.
  + [`fitSupportTimes`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/fitsupporttimes.html) Sorted unique observation times used as trajectory support.
  + [`mesoscaleConstraint`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/mesoscaleconstraint.html) Hard constraint applied to the fitted mesoscale streamfunction.
  + [`observedTrajectories`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/observedtrajectories.html) Observed drifter trajectory splines used for the fit.
  + [`psiKnotPoints`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/psiknotpoints.html) Mesoscale tensor-product knot vectors `{qKnot, rKnot, tKnot}`.
  + [`psiS`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/psis.html) Mesoscale spline degrees `[Sq Sr St]`.
  + [`representativeTimes`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/representativetimes.html) Representative pooled times from the stride-rule fast basis.
  + [`streamfunctionSpline`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/streamfunctionspline.html) Fitted COM-frame mesoscale streamfunction spline.
+ Inspect decomposition trajectories
  + [`decomposition`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/decomposition.html) Per-drifter decomposition trajectories in fixed and centered frames.
+ Apply fitted decomposition
  + [`decomposeTrajectories`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/decomposetrajectories.html) Apply the fitted decomposition to a supplied drifter ensemble.
+ Evaluate fitted mesoscale
  + [`psiMesoscale`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/psimesoscale.html) Evaluate the fitted mesoscale streamfunction.
  + [`uMesoscale`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/umesoscale.html) Evaluate the mesoscale x-velocity $$-\psi_{\tilde{y}}$$.
  + [`vMesoscale`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/vmesoscale.html) Evaluate the mesoscale y-velocity $$\psi_{\tilde{x}}$$.
+ Evaluate fitted diagnostics
  + [`centeredCoordinates`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/centeredcoordinates.html) Convert fixed-frame coordinates to COM-frame coordinates.
  + [`sigma_n`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/sigma_n.html) Evaluate the normal strain field $$\sigma_n = -2\psi_{\tilde{x}\tilde{y}}$$.
  + [`sigma_s`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/sigma_s.html) Evaluate the shear strain field $$\sigma_s = \psi_{\tilde{x}\tilde{x}} - \psi_{\tilde{y}\tilde{y}}$$.
  + [`uBackground`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/ubackground.html) Evaluate the fitted background x-velocity.
  + [`vBackground`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/vbackground.html) Evaluate the fitted background y-velocity.
  + [`zeta`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/zeta.html) Evaluate the relative-vorticity field $$\zeta = \psi_{\tilde{x}\tilde{x}} + \psi_{\tilde{y}\tilde{y}}$$.
+ Visualize strain angle
  + [`visualPrincipalStrainAngle`](/advection-diffusion-models/classes/estimators/gridded-streamfunction/griddedstreamfunction/visualprincipalstrainangle.html) Evaluate a jump-free principal strain angle for visualization.


---