---
layout: default
title: mesoscaleDegreesOfFreedom
parent: GriddedStreamfunction
grand_parent: Classes
nav_order: 12
mathjax: true
---

#  mesoscaleDegreesOfFreedom

Identifiable mesoscale degrees of freedom after gauge and constraints.


---

## Description
Real valued property with no dimensions and no units.

## Discussion

  `mesoscaleDegreesOfFreedom` counts the solved mesoscale spline
  degrees of freedom after removing the additive streamfunction
  gauge and applying the selected hard mesoscale constraint. This
  scalar is determined by the resolved mesoscale basis and the
  chosen structural constraint, not by the fitted coefficient
  values themselves.
