"""
This is in a seperate file to allow easy rerunning.
"""

using Gadfly: plot, layer, cm, Gadfly, Theme, Guide, Geom, Col, mm, style, Scale, PNG, SVG, SVGJS, Coord, draw
using Serialization, Cairo, DataFrames

# Gadfly defaults
Gadfly.set_default_plot_size(20cm, 20cm)
Gadfly.push_theme(Theme(background_color="white"))

include("utils.jl")

function populationCurves(result::Result, parameter::Parameters, plotFrame::Union{DataFrame, Nothing}=nothing)
  # Cell count curve
  counts = result.counts
  h = plot(counts, x=:time, y=:count, color=:genotype, Geom.line, Guide.title("$(parameter)"))
  if plotFrame ≠ nothing
    append!(plotFrame, DataFrame(plots=h, plotType="Population"))
  end
  return h
end

function proteinHistograms(result::Result, parameter::Parameters, plotFrame::Union{DataFrame, Nothing}=nothing)
  # Protein level histograms by genotype, time and protein
  levels = result.proteinLevels
  levels = levels[levels.level .> 0.0, :]
  genotypes = unique(levels.genotype)
  times = unique(levels.time)
  proteins = unique(levels.protein)
  for t in times
    for g in genotypes
      for p in proteins
        if size(levels[levels.protein .== p .&& levels.time .== t .&& levels.genotype .== g, :], 1) == 0
          levels = append!(levels, DataFrame(time=t, level=1.0, protein=p, genotype=g))
        end
      end
    end
  end

  layers = [
    layer(levels[levels.genotype .== gt, :], x=:level, Geom.histogram)
    for gt in genotypes
  ]
  h = plot(levels, 
    xgroup=:protein, 
    ygroup=:time,
    x=:level, 
    color=:genotype, 
    Geom.subplot_grid(layers..., Coord.cartesian(xmin=-1, xmax=1)),
    Guide.title("$(parameter)"),
    Scale.x_log10,
    )
  if plotFrame ≠ nothing
    append!(plotFrame, DataFrame(plots=h, plotType="Protein histograms"))
  end
  return h
end

function deathTimeHistograms(result::Result, parameter::Parameters, plotFrame::Union{DataFrame, Nothing}=nothing)
  # Death time histograms  
  deathTimes = result.deathTimes
  h = plot(x=deathTimes, color=:genotype, Geom.histogram(), Guide.xlabel("Age (hours)"), Guide.title("$(parameter)"))
  if plotFrame ≠ nothing
    append!(plotFrame, DataFrame(plots=h, plotType="Death time histograms"))
  end
  return h
end

if !isdefined(Main, :results)
  results = deserialize("results.dat");
end

plots = DataFrame()
for (parameter, result) in results
  local h
  local levels

  if parameter != ConcreteParameters(2.0, 0.3, "low BCLxL")
    continue
  end

  # h = populationCurves(result, parameter, plots)
  # display(h)
  # h |> PNG("/Users/thomas.e/Desktop/population $(parameter).png", 15cm, 15cm)

  try
    h = proteinHistograms(result, parameter, plots)
    h |> png
    # h |> SVG("/Users/thomas.e/Desktop/protein level $(parameter).svg", 40cm, 25cm)
  catch e
    @error "Fail" exception=(e, catch_backtrace())
  end
end


open("plots.dat", "w") do io
  serialize(io, plots)
end
