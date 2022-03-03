"""
This is in a seperate file to allow easy rerunning.
"""

for (parameter, result) in results
  # Cell count curve
  counts = result.counts
  h = plot(counts, x=:time, y=:count, Geom.line, Guide.title("$(parameter)"))
  display(h)

  # levels = result.proteinLevels
  # for protein in vcat(proteins, ["ensemble"])
  #   local ensembleLevel = levels[levels.name .== protein, :]
  #   local histoLayer = layer(ensembleLevel, x=:level, colour=:time, Geom.histogram(density=true))
  #   local densityLayer = layer(ensembleLevel, x=:level, Geom.density(bandwidth=0.01), style(line_width=0.8mm), color=:time)
  #   layers = [histoLayer, densityLayer]
  #   if protein == "ensemble"
  #     vline = layer(xintercept=[initialThreshold], Geom.vline(size=[1mm], color=["black"]))
  #     push!(layers, vline)
  #   end
  #   local h = Gadfly.plot(layers..., Guide.xlabel("$(protein) level"), Scale.x_log10)
  #   display(h)
  # end
end