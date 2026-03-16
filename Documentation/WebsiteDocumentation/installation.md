---
layout: default
title: Installation
nav_order: 2
description: Installation instructions"
permalink: /installation
---

## Installation

You may install this as part of [OceanKit](https://github.com/JeffreyEarly/OceanKit), or directly from [the github repo](https://github.com/JeffreyEarly/advection-diffusion-models).

### Package Manger Installation

Clone OceanKit
```
git clone https://github.com/JeffreyEarly/OceanKit.git
```
from the command-line. Within Matlab, add this folder as an MPM repository,
```matlab
mpmAddRepository("OceanKit","path/to/folder/OceanKit")
```
and then install the advection-diffusion-models package,
```matlab
mpminstall("advection-diffusion-models")
```
which will install the appropriate dependencies.

### Direct Installation

Installing the advection-diffusion-models package directly allows you to work directly from the github repo. You probably want to install the OceanKit mpm repository first, so that that the dependencies can be resolved. After that, you just need to clone the repo

```
https://github.com/JeffreyEarly/advection-diffusion-models
```
and then install
```matlab
mpminstall("local/path/to/advection-diffusion-models", Authoring=true);
```
with authoring enabled.
