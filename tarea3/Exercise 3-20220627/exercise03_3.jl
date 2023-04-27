# Setup and install packages
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


# Load packages
using Gridap
using GridapGmsh

using PrettyTables: pretty_table, ft_printf


# the analytical solution and its gradient
u(x) = 65 / 32 - (x[1]^2 + x[2]^2)  / 4 + 31 / (32 * log(4)) * log(x[1]^2 + x[2]^2)
function ∇u(x)
  factor = -1 / 2 + 31 / (16 * (x[1]^2 + x[2]^2) * log(4))
  return VectorValue(factor * x[1], factor * x[2])
end

# Inform Gridap.jl about the gradient of `u`
# see https://gridap.github.io/Tutorials/stable/pages/t002_validation/#Manufactured-solution-1
Gridap.∇(::typeof(u)) = ∇u

# the right-hand side
f(x) = 1


"""
    solve_poisson(mesh_file)

Solve the given Poisson equation on the mesh defined by `mesh_file`.
Returns the discrete L2 and H1 errors.
"""
function solve_poisson(mesh_file; order=1)
  # Load the mesh
  @info "Loading mesh file $mesh_file"
  model = GmshDiscreteModel(mesh_file)
  euler_characteristic = num_vertices(model) - num_edges(model) + num_cells(model)
  @info "num_vertices(model) - num_edges(model) + num_cells(model) = $euler_characteristic"

  # Set up the FE spaces
  reference_fe = ReferenceFE(lagrangian, Float64, order)
  V0 = TestFESpace(model, reference_fe; conformity=:H1,
                   dirichlet_tags=["inner_boundary", "outer_boundary"])
  U = TrialFESpace(V0, u)

  # Set up the triangulation and quadrature
  degree = 2 * order # to integrate the mass (and stiffness) matrix exactly
  Ω = Triangulation(model)
  dΩ = Measure(Ω, degree)

  # Set up the weak form
  a(u,v) = ∫( ∇(v) ⋅ ∇(u) ) * dΩ
  b(v) = ∫( v*f ) * dΩ
  op = AffineFEOperator(a, b, U, V0)

  # Solve the discrete problem
  uh = solve(op)

  # Compute the discrete errors
  e = u - uh
  error_l2 = sqrt(sum( ∫( e*e )*dΩ ))
  error_h1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩ ))

  # Write data to a file for manual inspection
  filename = first(splitext(mesh_file)) * "_solution.vtu"
  writevtk(Ω, filename, cellfields=["error" => e, "uh" => uh, "u" => u])

  return error_l2, error_h1
end


"""
    exercise03_3()

Perform the task of exercise 3.3 and print numerical errors as well as the
experimental order of convergence (EOC) to the screen.
"""
function exercise03_3(; order=1)
  # find all mesh files in the current directory
  meshes = String[]
  for (root, dirs, files) in walkdir(@__DIR__)
    for file in files
      if endswith(file, ".msh")
        push!(meshes, joinpath(root, file))
      end
    end
  end

  # setup arrays to store the results
  errors_l2 = zeros(length(meshes))
  errors_h1 = zeros(length(meshes))
  eoc_l2 = zeros(length(meshes))
  eoc_h1 = zeros(length(meshes))

  # compute errors in the maximum norm
  for (idx, mesh_file) in enumerate(meshes)
    err_l2, err_h1 = solve_poisson(mesh_file; order)
    errors_l2[idx] = err_l2
    errors_h1[idx] = err_h1
  end

  # compute experimental order of convergence for meshes with a refinement ratio of 2
  eoc_l2[1] = eoc_h1[1] = NaN # no EOC defined for the first grid
  for idx in 2:length(meshes)
    eoc_l2[idx] = -( log(errors_l2[idx] / errors_l2[idx - 1]) / log(2) )
    eoc_h1[idx] = -( log(errors_h1[idx] / errors_h1[idx - 1]) / log(2) )
  end

  # print results
  data = hcat(errors_l2, eoc_l2, errors_h1, eoc_h1)
  header = ["L2 error", "L2 EOC", "H1 error", "H1 EOC"]
  kwargs = (; header, formatters=(ft_printf("%.2e", [1, 3]), ft_printf("%.2f", [2, 4])))
  pretty_table(data; kwargs...)
  pretty_table(data; kwargs..., backend=Val(:latex))

  return nothing
end
