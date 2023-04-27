# Setup and install packages
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


# Load packages
using Gridap
using GridapGmsh

function save_visualization(mesh_file)
  @info "Loading mesh file $mesh_file"

  # Load the mesh
  model = GmshDiscreteModel(mesh_file)
  euler_characteristic = num_vertices(model) - num_edges(model) + num_cells(model)
  @info "num_vertices(model) - num_edges(model) + num_cells(model) = $euler_characteristic"

  writevtk(model, "model")
end
