"Model a distribution of time courses."

using Distributions, DataFrames, Fontconfig, Cairo
import Base.show, Base.append!

abstract type Parameters end

abstract type TimeCourseParameters end
function (f::TimeCourseParameters)(t::Float64) end

struct GammaTimeCourseParms <: TimeCourseParameters
  A::Float64
  α::Float64
  β::Float64
end
function GammaTimeCourseParms(d::Distribution, A::Float64, α::Float64, β::Float64)
  GammaTimeCourseParms(rand(d)*A, α, β)
end
function (f::GammaTimeCourseParms)(time::Float64)
  f.A * pdf.(Gamma(f.α, f.β), time)
end

struct PiecewiseLinear <: TimeCourseParameters
  A::Float64
  tMax::Float64
end
function PiecewiseLinear(d::Distribution, tMax::Float64)
  PiecewiseLinear(rand(d), tMax)
end
function (f::PiecewiseLinear)(time::Float64)
  tp = 72.0 # turning point
  if time < tp
    r = time/tp
  elseif time < f.tMax
    r = 1.0 - (time-tp)/(f.tMax-tp)
  else
    r = 0.0
  end

  return f.A * r
end

struct Result
  counts::DataFrame
  proteinLevels::DataFrame
  deathTimes::Vector{Float64}
end
Result() = Result(DataFrame(), DataFrame(), Float64[])

function Base.append!(r1::Result, r2::Result)
  append!(r1.counts, r2.counts)
  append!(r1.proteinLevels, r2.proteinLevels)
  append!(r1.deathTimes, r2.deathTimes)
  return r1
end

struct ConcreteParameters <: Parameters
  threshold::Float64
  gstd::Float64
  comment::String
end
ConcreteParameters(threshold::Float64, gstd::Float64) = ConcreteParameters(threshold, gstd, "")

show(io::IO, ::Parameters) = print(io, "no desciption")
function show(io::IO, parms::ConcreteParameters) 
  if parms.comment ≠ ""
    extra = " ($(parms.comment))"
  else
    extra = ""
  end
  print(io, "threshold=$(parms.threshold) gstd=$(parms.gstd)$extra")
end


# # Render to PNG instead of SVG
# # https://discourse.julialang.org/t/why-is-julias-graphics-system-so-slow/68750/33
# struct PNGPlot
#     p::Gadfly.Plot
# end
# Base.show(io::IO, ::MIME"image/png", pp::PNGPlot) = draw(PNG(io), pp.p)

# png(p) = PNGPlot(p)

#plot(y=[1,2,3]) |> png
