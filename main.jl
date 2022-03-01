using Cyton

println(rpad(lpad(" start ", 30, "-"), 55, "-"))

include("src/models/stretched.jl")

model = createModel(1000, stretchedCellFactory)

run(model, 500.0)
