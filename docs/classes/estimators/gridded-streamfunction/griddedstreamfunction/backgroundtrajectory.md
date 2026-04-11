---
layout: default
title: backgroundTrajectory
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  backgroundTrajectory

Fitted common background trajectory.


---

## Discussion

  `backgroundTrajectory.x(t)` evaluates $$x^{\mathrm{bg}}(t)$$
  and `backgroundTrajectory.y(t)` evaluates
  $$y^{\mathrm{bg}}(t)$$, with
  $$x^{\mathrm{bg}}(t_0)=0$$ and $$y^{\mathrm{bg}}(t_0)=0$$ at
  the global fit start time. The recovered background velocity is
  obtained from `backgroundTrajectory.u(t)` and
  `backgroundTrajectory.v(t)`.


