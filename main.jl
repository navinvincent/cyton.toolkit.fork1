using Agents, InteractiveDynamics, StatsPlots, DataFrames
plotlyjs()

using Cyton: CellAgent, step, VoidSpace

# include("src/models/michelle.jl")
include("src/models/simple.jl")

model_time(model::AgentBasedModel) = model.properties[:step_cnt] * model.properties[:Δt]

function init(; N=100)
  space = VoidSpace()
  scheduler = Schedulers.fastest
  properties = Dict(:step_cnt => 0, :Δt => 1.0)
  model = AgentBasedModel(CellAgent, space; properties, scheduler)

  for id in 1:N
    cell = cellFactory()
    agent = CellAgent(id, cell)
    add_agent_pos!(agent, model)
  end

  return model
end

# print("Time to initialise:")
# @time 
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
  counts = DataFrame(time=Int[], total=[], predivision=[], dividing=[], destiny=[], memory=[])
  for time in 1:200
    step!(model, stepper, model_stepper)

    cellAgents = values(model.agents)
    tm = model_time(model)
    local memoryCnt = 0
    local preDivisionCnt = 0
    local dividingCnt = 0
    local destinyCnt = 0
    for c in cellAgents
      ct = c.cell.differentiationAccumulator.cellType
      if ct == Memory
        memoryCnt += 1
        continue
      end
      if ct == Undivided
        preDivisionCnt += 1
        continue
      end
      if ct == Dividing
        dividingCnt += 1
        continue
      end
      if ct == Destiny
        destinyCnt += 1
        continue
      end
      error("Unkown cell type $(ct)")
    end
    local totalCnt = preDivisionCnt + dividingCnt + destinyCnt + memoryCnt
    push!(counts, (time, totalCnt, preDivisionCnt, dividingCnt, destinyCnt, memoryCnt))
  end
end


h = @df counts plot(:time, [:total :predivision :dividing :destiny :memory])
display(h)
println("Done at model time=$(model.properties[:step_cnt]*model.properties[:Δt])")

