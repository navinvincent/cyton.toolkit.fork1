using Agents, InteractiveDynamics, StatsPlots, DataFrames
plotlyjs()

using Cyton: createModel

function dashedLine(n::Int=20, c::Char='-')
  s = ""
  for _ in 1:n
    s *= c
  end
  return s
end
d = dashedLine()

println(d * " start " * d)

include("src/models/stretched.jl")

model = createModel(1000, stretchedCellFactory)

run(model, 500.0)
