---
layout: default
title: momentTensorEvolution
parent: LinearVelocityField
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  momentTensorEvolution

Evolve the second-moment tensor for the linear model.


---

## Declaration
```matlab
 [Mxx, Myy, Mxy] = momentTensorEvolution(self,Mxx0,Myy0,Mxy0,t,kappa)
```
## Parameters
+ `Mxx0`  initial $$M_{xx}$$ entry in m^2
+ `Myy0`  initial $$M_{yy}$$ entry in m^2
+ `Mxy0`  initial $$M_{xy}$$ entry in m^2
+ `t`  time vector in seconds
+ `kappa`  scalar diffusivity in $$m^2 s^-1$$

## Returns
+ `Mxx`  time-dependent $$M_{xx}$$ entry in m^2 with the same shape as `t`
+ `Myy`  time-dependent $$M_{yy}$$ entry in m^2 with the same shape as `t`
+ `Mxy`  time-dependent $$M_{xy}$$ entry in m^2 with the same shape as `t`

## Discussion

  The returned moments satisfy the matrix equation
  $$\dot{M} = AM + MA^\top + 2\kappa I$$ for the affine velocity gradient
  `A` defined by the model parameters.


