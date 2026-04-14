# OpenPDAC-13

[![SQAaaS badge](https://github.com/EOSC-synergy/SQAaaS/raw/master/badges/badges_150x116/badge_software_silver.png)](https://api.eu.badgr.io/public/assertions/fWEJSjZPRJ601jqAWQ9plg "SQAaaS silver badge achieved")

[![SQAaaS badge shields.io](https://img.shields.io/badge/sqaaas%20software-silver-lightgrey)](https://api.eu.badgr.io/public/assertions/fWEJSjZPRJ601jqAWQ9plg "SQAaaS silver badge achieved")

[![DOI](https://zenodo.org/badge/1016589416.svg)](https://doi.org/10.5281/zenodo.17054990)


OpenPDAC is an OpenFOAM module based on the module multiphaseEuler,
distributed with OpenFOAM.

With respect to the origianl module, in OpenPDAC the equations from
the kinetic theory for granular flows are modified to model multiple
dispersed solid phases.

In addition, a lagrangian library is included in the model (one-way
coupling with the gas-solid mixture).

The module also implement an initialization of the hydrostatic pressure
profile, which is needed for simulations on large domains. This allows
you to use boundary conditions which are appropriate for inflow/outflow.

Five test cases are provided:

- a 3D explosion simulation;
- a 2D explosion simulation on a flat topography;
- a 2D dilute flow over a wavy surface;
- a 2D fluidezed bed with two solid phases;
- a 2D impinging flow with two solid phases.

This version is based on the Ubuntu package:
openfoam13_20260212_amd64.deb

This code is not approved not endorsed by the OpenFOAM Foundation or
by ESI Ltd, the owner of OpenFOAM.
