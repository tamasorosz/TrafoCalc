# TrafoCalc

[![Build](https://github.com/tamasorosz/TrafoCalc/actions/workflows/ci.yml/badge.svg)](https://github.com/tamasorosz/TrafoCalc/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/tamasorosz/TrafoCalc/branch/master/graph/badge.svg?token=6SBI4COCOQ)](https://codecov.io/gh/tamasorosz/TrafoCalc)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/tamasorosz/TrafoCalc.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/tamasorosz/TrafoCalc/context:python)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/tamasorosz/TrafoCalc.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/tamasorosz/TrafoCalc/alerts/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?style=flat-square)](https://github.com/psf/black)

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

## Quickstart

The `\notes` library contains many sample designs in jupyter notebooks, which can be a good starting point to use this
package.

## References

[1] `Orosz, T., Borbély, B., Tamus, Z. Ádám  
Performance Comparison of Multi Design Method and Meta-Heuristic Methods for Optimal Preliminary Design of Core-Form Power Transformers , Periodica Polytechnica Electrical Engineering and Computer Science, 61(1), pp. 69–76, 2017. https://doi.org/10.3311/PPee.10207`
