using Distributions, Gadfly, DataFrames

d1 = Normal(0, 1)
d2 = Normal(1, 0.5)

x = rand(d1, 10000)
df = DataFrame(x=x, class="Class 1")

x = 0.5 .* rand(d2, 10000)

append!(df, DataFrame(x=x, class="Class 2"))

spike = -2.0 .* ones(500)
append!(df, DataFrame(x=spike, class="Class 2"))

h = plot(df[df.class .== "Class 1", :], 
  x=:x, 
  color=:class, 
  Geom.histogram(bincount=200), 
  Coord.cartesian(xmin=-4, xmax=4, ymin=0, ymax=1000)
  )
display(h)

h = plot(df,
  x=:x, 
  color=:class, 
  Geom.histogram(bincount=200), 
  Coord.cartesian(xmin=-4, xmax=4, ymin=0, ymax=1000)
  )
display(h)

layers = [layer(df[df.class .== c, :], x=:x, color=:class, Geom.histogram(bincount=200)) for c in ["Class 1", "Class 2"]]
h = plot(layers...,
  Coord.cartesian(xmin=-4, xmax=4, ymin=0, ymax=1000)
  )
display(h)



