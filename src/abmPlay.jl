using Revise

add_paths = ["probability", "cells", "space"]

for p in add_paths
  pp = "src/" * p
  if !(pp in LOAD_PATH)
    push!(LOAD_PATH, pp)
  end
end
