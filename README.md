# TrafoCalc

[![Build](https://github.com/tamasorosz/TrafoCalc/actions/workflows/ci.yml/badge.svg)](https://github.com/tamasorosz/TrafoCalc/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/tamasorosz/TrafoCalc/branch/master/graph/badge.svg?token=6SBI4COCOQ)](https://codecov.io/gh/tamasorosz/TrafoCalc)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?style=flat-square)](https://github.com/psf/black)
[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2Ftamasorosz%2FTrafoCalc.svg?type=shield)](https://app.fossa.com/projects/git%2Bgithub.com%2Ftamasorosz%2FTrafoCalc?ref=badge_shield)

## Goal of the Project

Design and Optimization of a power transformer needs to handle and calculate many required technical parameters
simoultanously, while the econoimically or technically optimal design was searched. This python package contains the
implementations of many analytical functions, which can be used to calculate the geometrical, electrical, mechanical or
thermal properties of three-phased power transformers.

The tool uses the implementation of a simple two-winding model calculations for the calculation of the transformer
parameters. Beside the analytical formulas, this project contains a 2D axisymmetric, parametric FEM simulation for
calculation of the short circuit impedance of a power transformer. The proposed methodology can be useful to understand
the principles of transformer design, or transformer professionals. Because, all of the proposed calculations validated
by measurements.

## Installation
The project uses python3.8.10 and poetry for package management.
You can install this python version by pyenv then you can install the selected package by simply the following command:

> poetry install

or the implemented transformer optimization packages accessible via the pip package:

> pip install trafocalc

## Quickstart

The `\notes` library contains many sample designs in jupyter notebooks, which can be a good starting point to use this
package.

## Publications and Example Calculations

You can find the detailed description of the transformer optimization functions and example calculation in the
following paper:

```
Orosz, Tamás
FEM-Based Power Transformer Model for Superconducting and Conventional Power Transformer Optimization"
Energies 2022, 15(17), 6177; https://doi.org/10.3390/en15176177 
```

## Additional References

[2] `Orosz, T., Borbély, B., Tamus, Z. Ádám  
Performance Comparison of Multi Design Method and Meta-Heuristic Methods for Optimal Preliminary Design of Core-Form Power Transformers , Periodica Polytechnica Electrical Engineering and Computer Science, 61(1), pp. 69–76, 2017. https://doi.org/10.3311/PPee.10207`


## License
[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2Ftamasorosz%2FTrafoCalc.svg?type=large)](https://app.fossa.com/projects/git%2Bgithub.com%2Ftamasorosz%2FTrafoCalc?ref=badge_large)