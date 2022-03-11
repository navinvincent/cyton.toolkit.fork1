# Parameters for the distributions protein level time courses

using Gadfly: plot, layer, cm, Gadfly, Theme, Guide, Geom, Col, mm, style, Scale, PNG
using DataFrames, Colors, Cairo
using Base.Threads, Serialization
import Base.show

# Gadfly defaults
Gadfly.set_default_plot_size(20cm, 20cm)
Gadfly.push_theme(Theme(background_color="white"))

include("timecourses.jl")
TimeCourseParms(gstd::Float64) = GammaTimeCourseParms(LogNormal(0, gstd), 114.2, 4.6, 24.1)
# TimeCourseParms(gstd::Float64) = PiecewiseLinear(LogNormal(0, gstd), 235.0)

include("MiCell.jl")

proteins = ["BCL2", "BCLxL", "MCL1", "BIM"]
initialWeights = Dict{String, Float64}([
  "BCL2"  =>  5,
  "BCLxL" =>  5,
  "MCL1"  =>  5,
  "BIM"   => -5,
])
initialThreshold = 6.0
initialGstd = 0.3

struct ConcreteParameters <: Parameters
  threshold::Float64
  gstd::Float64
end
show(io::IO, ::Parameters) = print(io, "no desciption")
show(io::IO, parms::ConcreteParameters) = print(io, "threshold=$(parms.threshold) gstd=$(parms.gstd)")

thresholds = [8.0]#1.0:8.0
gstds      = [0.5]#0.1:0.1:0.5
parameters = ConcreteParameters[]
for threshold in thresholds
  for gstd in gstds
    push!(parameters, ConcreteParameters(threshold, gstd))
  end
end


snapshotTimes = [72.0, 96.0, 120.0, 144.0, 168.0]
ensembleSnapshots = DataFrame(time=Vector{String}(), ensemble=Vector{Float64}(), threshold=Vector{Float64}(), gstd=Vector{Float64}())
deathSnapshots = DataFrame(ages=Vector{Float64}(), threshold=Vector{Float64}(), gstd=Vector{Float64}())
lk = ReentrantLock()

@threads for parameter in parameters

  N = 10000
  timers = Vector{Union{ThresholdDeath, Nothing}}(undef, N)
  deathTimes = Vector{Float64}()
  sizehint!(deathTimes, N)

  for i in 1:N
    d = ThresholdDeath(parameter.threshold, initialWeights, Γ_timeCourseParms(parameter.gstd))
    timers[i] = d
  end

  runDuration = 200.0
  Δt = 0.1

  for time in 0.0:Δt:runDuration
    for i in 1:N
      d = timers[i]
      if d === nothing
        continue
      end
      if !d.isDead && shouldDie(d, time)
        push!(deathTimes, time)
        # timers[i] = nothing
      end
    end

    lock(lk) do
      append!(deathSnapshots, DataFrame(ages=deathTimes, threshold=parameter.threshold, gstd=parameter.gstd))
    end

    if time in snapshotTimes
      ensembles = [ensemble(d, time) for d in timers if d !== nothing]
      snapshot = DataFrame(time=string(time), ensemble=ensembles, threshold=parameter.threshold, gstd=parameter.gstd)
      lock(lk) do
        append!(ensembleSnapshots, snapshot)
      end
    end

  end

  # hh(tm) = layer(snapshots[snapshots.time.==string(tm), :], x=:ensemble, colour=:time, Geom.histogram(density=false))
  # histos = [hh(tm) for tm in snapshotTimes]
  # ll(tm) = layer(snapshots[snapshots.time.==string(tm), :], x=:ensemble, Geom.density(bandwidth=0.01), style(line_width=0.8mm), color=[colorant"black"])
  # densities = [ll(tm) for tm in snapshotTimes]
  # h = Gadfly.plot(
  #   histos...,
  #   densities...,
  #   Guide.xlabel("Ensemble level"), 
  #   # Scale.x_log10
  #   )
  # display(h)

  # layers = [
  #   layer(snapshots, ygroup=:time, x=:ensemble, Geom.histogram(density=false)),
  #   layer(snapshots, ygroup=:time, x=:ensemble, Geom.density(), style(line_width=0.8mm), color=[colorant"black"]),
  # ]
  # h = plot(
  #   Geom.subplot_grid(layers...),
  #   Guide.xlabel("Ensemble level"), 
  #   Scale.x_log10
  #   )
  # display(h)

end

open("ensembleSnapshots.dat", "w") do io
  serialize(io, ensembleSnapshots)
end
open("deathSnapshots.dat", "w") do io
  serialize(io, deathSnapshots)
end



for parameter in parameters
  # flabels = DataFrame(labels=["$(t)h" for t in snapshotTimes], time=snapshotTimes, x=100, y=0.5)
  th = parameter.threshold
  gstd = parameter.gstd
  snapshot = ensembleSnapshots[ensembleSnapshots.threshold .== th .&& ensembleSnapshots.gstd .== gstd, :]
  layers = [
    # layer(flabels,   ygroup=:time, x=:x, y=:y, label=:labels, Geom.label(position=:centered)),
    layer(snapshot, ygroup=:time, x=:ensemble, Geom.histogram(density=true)),
    layer(snapshot, ygroup=:time, x=:ensemble, Geom.density(), style(line_width=0.8mm), color=[colorant"black"]),
    ]
  h = plot(
    Geom.subplot_grid(layers...),
    Guide.xlabel("Ensemble levels for $(parameter)"), 
    Scale.x_log10
    )
  display(h)
  h |> PNG("/Users/thomas.e/Desktop/gamma time course/ensemble histograms $(parameter).png", 15cm, 15cm)

  deathTimes = deathSnapshots[deathSnapshots.threshold .== th .&& deathSnapshots.gstd .== gstd, :ages]
  h = plot(x=deathTimes, Geom.histogram(), Guide.xlabel("Age (hours)"), Guide.title("$(parameter)"))
  display(h)
  h |> PNG("/Users/thomas.e/Desktop/gamma time course/death times $(parameter).png", 15cm, 15cm)

end