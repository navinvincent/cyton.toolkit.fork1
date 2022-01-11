module Cells
using Agents, DistributionParms

export Cell

mutable struct Cell <: AbstractAgent
  id::Int
  pos::NTuple{2, Int}
  d::DistributionParmSet
  λ::Float64
  start::Float64
  deaths::Int
  function Cell(id::Int, pos::NTuple{2, Int}, d::DistributionParmSet)
    return new(id, pos, d, draw(d), 0)
  end
end

end