import Agents: 
  AgentBasedModel, 
  AbstractAgent, 
  AbstractSpace, 
  random_position, 
  move_agent!,
  add_agent_to_space!, 
  remove_agent_from_space!, 
  nearby_ids

struct VoidSpace <: AbstractSpace end

random_position(model::AgentBasedModel{<:VoidSpace}) = (0, 0)
move_agent!(agent::AbstractAgent, pos, model::AgentBasedModel{<:VoidSpace}) = nothing
add_agent_to_space!(agent::AbstractAgent, model::AgentBasedModel{<:VoidSpace}) = nothing 
remove_agent_from_space!(agent::AbstractAgent, model::AgentBasedModel{<:VoidSpace}) = nothing
nearby_ids(position, model::AgentBasedModel{<:VoidSpace}, r) = Vector{NTuple{2, Int}}()

# add_agent_single!(agent::A, model::AgentBasedModel{<:VoidSpace}) where {A<:AbstractAgent} = agent

