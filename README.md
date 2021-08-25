# ContinuousTimeMarkov.jl

![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
[![build](https://github.com/tpapp/ContinuousTimeMarkov.jl/workflows/CI/badge.svg)](https://github.com/tpapp/ContinuousTimeMarkov.jl/actions?query=workflow%3ACI)
[![codecov.io](http://codecov.io/github/tpapp/ContinuousTimeMarkov.jl/coverage.svg?branch=master)](http://codecov.io/github/tpapp/ContinuousTimeMarkov.jl?branch=master)

## Overview

A small package for working with continuous-time Markov transition rate matrices.

Currently exposes the following functionality:

- wrappers `TransitionRateMatrix` and `TransitionRateMatrix!` for creating transition rate matrices (square, nonnegative elements, rows sum to one),

- `steady_state` for calculating steady state distributions.

## References

- Meyer, Carl D. "Sensitivity of the stationary distribution of a Markov chain." *SIAM Journal on Matrix Analysis and Applications* 15.3 (1994): 715-728.
