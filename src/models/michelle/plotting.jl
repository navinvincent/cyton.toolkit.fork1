"""
This is in a seperate file to allow easy rerunning.
"""

using Gadfly: plot, layer, cm, Gadfly, Theme, Guide, Geom, Col, mm, style, Scale, PNG

if !isdefined(Main, :results)
  results = deserialize("results.dat");
end

for (parameter, result) in results
  # Cell count curve
  counts = result.counts
  h = plot(counts, x=:time, y=:count, Geom.line, Guide.title("$(parameter)"))
  display(h)
  # h |> PNG("/Users/thomas.e/Desktop/population $(parameter).png", 15cm, 15cm)

  # # Death time histograms
  # deathTimes = result.deathTimes
  # h = plot(x=deathTimes, Geom.histogram(), Guide.xlabel("Age (hours)"), Guide.title("$(parameter)"))
  # display(h)
  # # h |> PNG("/Users/thomas.e/Desktop/death histogram $(parameter).png", 15cm, 15cm)

  # # Protein level histograms
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