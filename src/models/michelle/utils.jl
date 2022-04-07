"Model a distribution of time courses."

using Distributions, DataFrames, Fontconfig, Cairo
import Base.show, Base.append!

abstract type Parameters end
abstract type JobResult end

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
  bclxlWeight::Float64
  inhibitionFactor::Float64
  cellType::CellType
  comment::String
end
ConcreteParameters(threshold::Float64, gstd::Float64, weight::Float64, inhibitionFactor::Float64, cellType::CellType) = ConcreteParameters(threshold, gstd, weight, inhibitionFactor, cellType, "")

struct ParameterKey <: Parameters
  threshold::Float64
  gstd::Float64
  bclxlWeight::Float64
  inhibitionFactor::Float64
end
ParameterKey(cp::ConcreteParameters) = ParameterKey(cp.threshold, cp.gstd, cp.bclxlWeight, cp.inhibitionFactor)

function parameterKey(::Parameters) end
parameterKey(cp::ConcreteParameters) = ParameterKey(cp)

show(io::IO, ::Parameters) = print(io, "no desciption")
function show(io::IO, parms::ParameterKey) 
  print(io, "threshold=$(parms.threshold) gstd=$(parms.gstd) BCLxL weight=$(parms.bclxlWeight) IF=$(parms.inhibitionFactor)")
end
function show(io::IO, parms::ConcreteParameters) 
  if parms.comment == ""
    extra = ""
  else
    extra = " ($(parms.comment))"
  end
  print(io, "threshold=$(parms.threshold) gstd=$(parms.gstd) BCLxL weight=$(parms.bclxlWeight) IF=$(parms.inhibitionFactor), cell type=$(parms.cellType)$extra")
end
