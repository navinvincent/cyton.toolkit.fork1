"""
This is in a seperate file to allow easy rerunning.
"""

using Gadfly: plot, layer, cm, Gadfly, Theme, Guide, Geom, Col, mm, style, Scale, PNG, SVG, SVGJS, Coord
using Serialization, Cairo, DataFrames, Base.Threads

# Gadfly defaults
Gadfly.set_default_plot_size(20cm, 20cm)
Gadfly.push_theme(Theme(background_color="white", alphas=[0.5]))

include("utils.jl")

function populationCurves(counts::DataFrame, y::Symbol, title::String)
  # Cell count curve
  return plot(counts, 
                x=:time, 
                y=y, 
                color=:genotype, 
                Geom.line, 
                style(line_width=1mm),
                Guide.title(title))
end

function proteinHistograms2(levels::DataFrame, title::String)
  levels = levels[levels.level .> 0.0, :]
  h = plot(levels,
    x=:level,
    color=:genotype,
    Geom.histogram(position=:dodge),
    Coord.cartesian(xmin=-1, xmax=1),
    Guide.title(title),
    Scale.x_log10,
    Scale.x_log10(minvalue=0.1, maxvalue=10),
    )
  return h
end
  
function proteinHistograms(levels::DataFrame, title::String)
  # Protein level histograms by genotype, time and protein
  l = levels
  levels = levels[l.level .> 0.0 .&& l.time.>=140 .&& l.time.<=180 .&& l.protein.!="ensemble", :]
  if isempty(levels)
    return nothing
  end
  h = plot(levels,
    xgroup=:protein,
    ygroup=:time,
    x=:level,
    color=:genotype,
    Geom.subplot_grid(Geom.histogram(position=:dodge), Coord.cartesian(xmin=-1, xmax=0.5)),
    Guide.title(title),
    Scale.x_log10
    )
  return h
end

function deathTimeHistograms(result::Result, parameter::Parameters)
  # Death time histograms  
  deathTimes = result.deathTimes
  h = plot(x=deathTimes, color=:genotype, Geom.histogram(), Guide.xlabel("Age (hours)"), Guide.title("$(parameter)"))
  return h
end

if !isdefined(Main, :results)
  results = deserialize("results.dat");
end

function bigPlot(parameter, result)
  counts = result.counts
  levels = result.proteinLevels
  fileCheck = false

  fn = "outputs/population $(parameter).png"
  if !fileCheck || !isfile(fn)
    h = populationCurves(counts, :count, "$(parameter)")
    # display(h)
    h |> PNG(fn, 15cm, 15cm)
  end
  
  fn = "outputs/cohort $(parameter).png"
  if !fileCheck || !isfile(fn)
    h = populationCurves(counts, :cohort, "$(parameter)")
    # display(h)
    h |> PNG(fn, 15cm, 15cm)
  end
  
  # fn = "outputs/protein level $(parameter).png"
  # if !fileCheck || !isfile(fn)
  #   h = proteinHistograms(levels, "$(parameter)")
  #   if !isnothing(h)
  #     # display(h)
  #     h |> PNG(fn, 40cm, 25cm)
  #   end
  # end
  
  # times = unique(levels[!, :time])
  # local proteins = unique(levels[!, :protein])
  # for time in times
  #   for protein in proteins
  #     title = "$(parameter) time=$(time) protein=$(protein)"
  #     fn = "outputs/protein level $(title).png"
  #     if !fileCheck || !isfile(fn)
  #       lvl = levels[levels.time .== time .&& levels.protein .== protein, :]
  #       h = proteinHistograms2(lvl, title)
  #       # display(h)
  #       h |> PNG(fn, 15cm, 15cm)
  #     end
  #   end
  # end

end

struct PlotParameters <: Parameters
  actualParameters::ParameterKey
  result::Result
end
