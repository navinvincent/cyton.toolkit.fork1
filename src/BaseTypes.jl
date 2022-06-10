using Agents

export Time, Duration, CytonAgent
export CytonAgent, GlobalAgent, EnvironmentalAgent
export CytonEvent, GlobalEvent
export AgentImpl

#----------------------- Time based types -----------------------
"""
Time

A type alias for time-like quantities. Used to make function signatures more readable
"""
const Time = Float64

"""
Duration

A type alias for duration-like quantities. Used to make function signatures more readable
"""
const Duration = Float64
#----------------------------------------------------------------


#---------------------- Basic agent types -----------------------
"""
CytonAgent

The base type in the API for a Cyton agent. This can be a cell of another
agent that interacts with the cell population.
"""
abstract type CytonAgent end

"""
GlobalAgent

A type to model agents that effect a large number or all cells. Events generated by global
agents are propogated to all cells and global agents 
"""
abstract type GlobalAgent <: CytonAgent end

"""
EnvironmentalAgent

A type to model things in the cellular environment like cytokines, etc 
"""
abstract type EnvironmentalAgent <: GlobalAgent end
#----------------------------------------------------------------


#---------------------- Basic event types -----------------------
"""
CytonEvent

The base event type
"""
abstract type CytonEvent end

"""
GlobalEvent

An event that all agents potentially want to listen for
"""
abstract type GlobalEvent <: CytonEvent end
#----------------------------------------------------------------

#-------- Types for wrapping the Agent framework ----------
# This is used internally to wrap the Agents.jl type and is not exposed through the API.
mutable struct AgentImpl <: AbstractAgent
  "These are required by the ABM framework"
  id::Int
  pos::NTuple{2, Int}
  agent::CytonAgent
end
AgentImpl(agent::CytonAgent) = AgentImpl(0, (0, 0), agent)
AgentImpl(id::Int, agent::CytonAgent) = AgentImpl(id, (0, 0), agent)
#----------------------------------------------------------
