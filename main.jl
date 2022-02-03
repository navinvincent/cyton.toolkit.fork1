using Agents, InteractiveDynamics, StatsPlots, DataFrames
plotlyjs()

using Cyton: CellAgent, step, VoidSpace

function dashedLine(n::Int=20, c::Char='-')
  s = ""
  for _ in 1:n
    s *= c
  end
  return s
end
d = dashedLine()

println(d * " start " * d)

# include("src/models/michelle.jl")
include("src/models/simple.jl")

model_time(model::AgentBasedModel) = model.properties[:step_cnt] * model.properties[:Δt]

function init(; N=7477)
  space = VoidSpace()
  scheduler = Schedulers.fastest
  properties = Dict(:step_cnt => 0, :Δt => 0.1)
  model = AgentBasedModel(CellAgent, space; properties, scheduler)

  for id in 1:N
    cell = cellFactory()
    agent = CellAgent(id, cell)
    add_agent_pos!(agent, model)
  end

  return model
end

model = init()

function stepper(agent::CellAgent, model::AgentBasedModel)
  Δt = model.properties[:Δt]
  time = model_time(model)
  step(agent, time, Δt, model)
end

function model_stepper(model)
  model.properties[:step_cnt] += 1
end

print("Time to run:")
@time begin
  counts = DataFrame(time=Float64[], 
  total=[], 
  gen0 = [],
  gen1 = [],
  gen2 = [],
  gen3 = [],
  gen4 = [],
  gen5 = [],
  gen6 = [],
  gen7 = [],
  gen8 = [],
  genOther = []
  )
  cellAgents = values(model.agents)
  Δt = model.properties[:Δt]
  for time in 1:Δt:150 
    step!(model, stepper, model_stepper)

    tm = model_time(model)
    local genCnts = zeros(10)
    for c in cellAgents
      gen = c.cell.generation
      if gen <= 8
        genCnts[gen+1] += 1
      else
        genCnts[10] += 1
      end
    end
    push!(counts, (tm, length(cellAgents), genCnts...))
  end
end

h = @df counts plot(:time, [:total :gen0 :gen1 :gen2 :gen3 :gen4 :gen5 :gen6 :gen7 :gen8 :genOther])
display(h)

println("Done at model time=$(model.properties[:step_cnt]*model.properties[:Δt])")

