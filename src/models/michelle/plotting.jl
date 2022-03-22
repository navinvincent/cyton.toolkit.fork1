"""
This is in a seperate file to allow easy rerunning.
"""

using Gadfly: plot, layer, cm, Gadfly, Theme, Guide, Geom, Col, mm, style, Scale, PNG, SVG, SVGJS, Coord
using Serialization, Cairo, DataFrames, Base.Threads

# Gadfly defaults
Gadfly.set_default_plot_size(20cm, 20cm)
Gadfly.push_theme(Theme(background_color="white"))

include("utils.jl")

function populationCurves(result::Result, parameter::Parameters)
  # Cell count curve
  counts = result.counts
  h = plot(counts, x=:time, y=:count, color=:genotype, Geom.line, Guide.title("$(parameter)"))
  return h
end

function proteinHistograms(result::Result, parameter::Parameters)
  # Protein level histograms by genotype, time and protein
  levels = result.proteinLevels
  levels = levels[levels.level .> 0.0, :]
  h = plot(levels,
    xgroup=:protein,
    ygroup=:time,
    x=:level,
    color=:genotype,
    Geom.subplot_grid(Geom.histogram(position=:dodge), Coord.cartesian(xmin=-1, xmax=1)),
    Guide.title("$(parameter)"),
    Scale.x_log10(minvalue=0.1, maxvalue=100),
    )
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

lk = ReentrantLock()
plots = DataFrame()
function cb(h) 
  lock(lk) do
    append!(plots, DataFrame(plots=h, plotType="Protein histograms"))
  end
end

#@threads 
for (parameter, result) in collect(results)

  h = populationCurves(result, parameter)
  # display(h)
  cb(h)
  h |> PNG("/Users/thomas.e/Desktop/population $(parameter).png", 15cm, 15cm)

  h = proteinHistograms(result, parameter)
  # display(h)
  cb(h)
  h |> PNG("/Users/thomas.e/Desktop/protein level $(parameter).png", 40cm, 25cm)
end

open("plots.dat", "w") do io
  serialize(io, plots)
end
