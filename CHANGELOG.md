# Version History

## [2.1.0] - 2026-04-09
- updated the linear-velocity estimation path to the `Distributions` 2.0 named-argument API
- raised the package dependency floor to `Distributions ^2.0`
- migrated the remaining spline call sites to the `SplineCore 2.1` constructor and factory APIs, including the gridded-streamfunction workflows and legacy second-moment helper paths

## [2.0.0] - 2026-03-25
- updated the package for `SplineCore` 2.0 compatibility, including the internal spline call sites and package dependency metadata
- modernized the `Integrators` subsystem with updated constructors, lowerCamel public methods, reorganized class folders, refreshed examples, and expanded inline API documentation
- restored the standard OceanKit website scaffolding with generated class reference pages for the four integrator classes
