using Agents, InteractiveDynamics, Plots
plotlyjs()

using Cyton: CellAgent, step

include("src/cells/michelle.jl")

model_time(model::AgentBasedModel, dt::Float64) = model.properties["step_cnt"] * dt 

function init(; N=10)
  space = nothing
  scheduler = Schedulers.fastest
  properties = Dict("step_cnt" => 0)
  model = AgentBasedModel(CellAgent, space; properties, scheduler)

  for n in 1:N
    for m in 1:N
      cell = createMichelleCell()
      agent = CellAgent(cell)
      add_agent_pos!(agent, model)
    end
  end

  return model
end

print("Time to initialise:")
@time model = init()

dt = 0.1
function stepper(agent::CellAgent, model::AgentBasedModel)
  time = model_time(model, dt)
  step(agent, time, dt, model)
end

function model_stepper(model)
  model.properties["step_cnt"] += 1
end

print("Time to run:")
#@time 
step!(model, stepper, model_stepper, 200)

print("Time to plot: ")
@time begin
  cells = values(model.agents)
  tm = model_time(model, dt)
  lifetimes = collect(remaining(c, tm) for c in cells)
  lifetimes = sort(lifetimes)
  ds = 1 .- (1:1:length(lifetimes))/length(lifetimes)
  h = plot(lifetimes, ds)
  display(h);
end

println("Done at model time=$(tm)")

