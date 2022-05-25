# Cyton API
These are the types and methods that implement the Cyton framework.

See this guide for an overview: [cyton-dynamics.github.io/](https://cyton-dynamics.github.io/)


```@meta
CurrentModule = Cyton
```

```@docs
CellPopulation
```

```@docs
step(::FateTimer, ::Time, ::Duration)
```

```@docs
inherit
```

```@docs
modelTime(model::CellPopulation)
```

```@docs
modelTimeStep(model::CellPopulation)
```

```@docs
cellCount(model::CellPopulation)
```

```@docs
createPopulation(nCells::Int, 
  cellFactory::Function; 
  eventCallbacks::Vector{Function}=Function[])
```

```@docs
step(model::CellPopulation, stimulus::T) where T<:Stimulus
```

```@docs
step(model::CellPopulation, stimuli::Vector{T}=Vector{Stimulus}()) where T<:Stimulus
```

```@docs
step(agent::CellAgent, model::CellPopulation, stimuli::Vector{T}) where T<:Stimulus
```

```@docs

```