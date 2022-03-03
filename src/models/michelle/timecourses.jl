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
