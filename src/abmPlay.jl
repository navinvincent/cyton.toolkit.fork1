
add_paths = ["src/probability", "src/cells"]

for p in add_paths
  if !(p in LOAD_PATH)
    push!(LOAD_PATH, p)
  end
end
