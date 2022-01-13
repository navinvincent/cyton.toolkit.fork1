using Agents, InteractiveDynamics, Plots, Cyton
plotlyjs()

model_time(model::AgentBasedModel, dt::Float64) = model.properties["step_cnt"] * dt 

function init(; N=100, distributionParms=LogNormalParms(1.5, 0.3))
  space = nothing
  scheduler = Schedulers.fastest
  properties = Dict("step_cnt" => 0)
  model = AgentBasedModel(Cell, space; properties, scheduler)

  for n in 1:N
    for m in 1:N
      agent = Cell((n-1)*N+m, (n, m), distributionParms)
      add_agent_single!(agent, model)
    end
  end

  return model
end

print("Time to initialise:")
@time model = init()

dt = 0.1
function stepper(agent::Cell, model::AgentBasedModel)
  time = model_time(model, dt)
  step_cell(agent, time, model)
end

function model_stepper(model)
  model.properties["step_cnt"] += 1
end

print("Time to run:")
@time step!(model, stepper, model_stepper, 200)

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

