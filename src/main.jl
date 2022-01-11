include("abmPlay.jl")

using Agents, InteractiveDynamics, Plots, DistributionParms, Cells
plotlyjs()
Plots.PlotlyJSBackend()


model_time(model::AgentBasedModel, dt::Float64) = model.properties["step_cnt"] * dt 

function init(; N=1000, distributionParms=LogNormalParms(1.5, 0.3))
  space = GridSpace((N, N); periodic = false)
  scheduler = Schedulers.randomly
  properties = Dict("step_cnt" => 0)
  model = AgentBasedModel(Cell, space; properties, scheduler)

  for n in 1:N
    for m in 1:N
      agent = Cell((n-1)*N+m, (n, m), distributionParms)
      add_agent_pos!(agent, model)
    end
  end

  return model
end

@time model = init()

dt = 0.1
function stepper(agent::Cell, model::AgentBasedModel)
  time = model_time(model, dt)
  if time > agent.λ + agent.start
    agent.start = time
    agent.λ = draw(agent.d)
    agent.deaths += 1
  end
end

function model_stepper(model)
  model.properties["step_cnt"] += 1
end
println("Init done")

println("Start run")
@time step!(model, stepper, model_stepper, 1000)
println("Finish run")

# function cellcolour(agent::Cell)
#   N = 100
#   cmap = range(HSL(colorant"red"), stop=HSL(colorant"green"), length=N)
#   λ_max = agent.d.useful_max
#   λ = min(agent.λ, λ_max) * (N-1) / λ_max
#   return cmap[Int(floor(λ+1))]
# end

# fig, _ = abm_plot(model; ac=cellcolour)
# display(fig)

cells = model.agents
lifetimes = Array{Float64, 1}(undef, length(cells))
time = model_time(model, dt)
for i in eachindex(cells)
  cell = cells[i]
  lifetimes[i] = cell.λ - (time - cell.start)
end

println("trying to plot")
h = histogram(lifetimes)
display(h);

println("Done at model time=$(time)")

