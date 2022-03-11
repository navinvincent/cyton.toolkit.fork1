"Model a distribution of time courses."

using Distributions

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