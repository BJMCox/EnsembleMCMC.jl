# EnsembleMCMC.jl

[![CI](https://github.com/BJMCox/EnsembleMCMC.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/BJMCox/EnsembleMCMC.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/BJMCox/EnsembleMCMC.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/BJMCox/EnsembleMCMC.jl)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://bjmcox.github.io/EnsembleMCMC.jl/)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE.md)

Ensemble MCMC with Stretch, differential-evolution, and snooker moves.
Supports fixed move mixtures, threaded evaluation, and resumable sampling.

```julia
using EnsembleMCMC, Random, Random123

rng = Philox4x((42, 1))
initial = [randn(rng, 2) for _ in 1:24]
state = initialize(rng, x -> -sum(abs2, x) / 2, initial)
step!(state, 100)               # warmup, without collecting draws
draws = sample!(state, 1_000)
draws.positions                  # coordinates × walkers × sweeps
```

[Documentation](https://bjmcox.github.io/EnsembleMCMC.jl/) · [Apache 2.0 license](LICENSE.md)
