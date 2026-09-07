# EnsembleMCMC.jl

Ensemble MCMC with Stretch, differential-evolution, and snooker moves.
Supports fixed move mixtures, threaded evaluation, and resumable sampling.

```julia
using EnsembleMCMC, Random, Random123

rng = Philox4x((42, 1))
initial = [randn(rng, 2) for _ in 1:24]
state = initialize(rng, x -> -sum(abs2, x) / 2, initial)
draws = sample!(state, 1_000)  # coordinates × walkers × sweeps
```

[Documentation](docs/src/index.md) · [Apache 2.0 license](LICENSE.md)
