# API

```@meta
CurrentModule = EnsembleMCMC
```

## Sampling

```@docs
initialize
step!
sample!
current_state
snapshot
```

## Moves

```@docs
StretchMove
DEMove
DESnookerMove
MoveMixture
```

## Execution

```@docs
SerialExecutor
ThreadedExecutor
```
