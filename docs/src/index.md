# Cyton API
These are the types and methods that implement the Cyton framework.

See this guide for an overview: [cyton-dynamics.github.io/](https://cyton-dynamics.github.io/)


```@meta
CurrentModule = Cyton
```

## Types and constructors
```@docs
Time
```

```@docs
Duration
```

### Events
```@docs
CellEvent
```

```@docs
Division
```

```@docs
Death
```

### Cells and populations
```@docs
CellPopulation
```

```@docs
CellType
```

```@docs
GenericCell
```

```@docs
AbstractCell{T} where T <: CellType
```

```@docs
Cell
```

```@docs
Cell(birth::Time)
```

```@docs
Cell(birth::Time, cellType::T) where T <: CellType
```

```@docs
Cell(birth::Time, divisionCount::Int64)
```

### Stimuli
```@docs
Stimulus
```

## API methods and functions
These are the methods than need to be implemented to build a model.

```@docs
stimulate
```

```@docs
step(::FateTimer, ::Time, ::Duration)
```

```@docs
inherit
```

## Utility methods
```@docs
addTimer
```

```@docs
age
```

```@docs
addObserver
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
cohortCount(model::CellPopulation)
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

### Probability distributions
```@docs
DistributionParmSet
```

```@docs
PoissonParms
```

```@docs
FixedDistributionParms
```

```@docs
draw(d::DistributionParmSet)
```

```@docs
usefulMax
```

```@docs
plotPdf
```

```@docs
NormalParms
```

```@docs
LogNormalParms
```
