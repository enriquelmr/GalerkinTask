# Setup and install packages
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


# Load packages from the standard library
using Printf
using SparseArrays

# Load a package from the general registry
using PrettyTables: pretty_table, ft_printf


"""
    setup(xmin, xmax, nx)

Setup the grid `x` of interior nodes and the stiffness matrix `A` for the
Laplace equation with homogeneous boundary conditions on a uniform grid in
the 1D domain `[xmin, xmax]` with `nx` interior nodes.

Returns `x, A`.
"""
function setup(xmin, xmax, nx)
  # Create a grid with two additional points and discard the endpoints
  x = range(xmin, xmax, nx+2)[(begin+1):(end-1)]

  h = step(x)
  bidiagonal = (-1 / h^2) * ones(nx - 1)
  diagonal = (2  / h^2) * ones(nx)
  A = spdiagm(-1 => bidiagonal,
              0  => diagonal,
              +1 => bidiagonal)

  return x, A
end


"""
    solve_exercise01_3a(nx)

Compute a numerical solution `uh` of exercise 1.3a with `nx` interior nodes of
the grid `x`.

Returns `x, uh`.
"""
function solve_exercise01_3a(nx)
  xmin = 0.0
  xmax = 1.0
  x, A = setup(xmin, xmax, nx)

  f = @. π^2 * sinpi(x)

  uh = A \ f
  return x, uh
end


"""
    solve_exercise01_3b(nx)

Compute a numerical solution `uh` of exercise 1.3b with `nx` interior nodes of
the grid `x`.

Returns `x, uh`.
"""
function solve_exercise01_3b(nx)
  xmin = 0.0
  xmax = 1.0
  x, A = setup(xmin, xmax, nx)

  # We subtract `w = 1 + (exp(1) - 1) * x` to get a problem with homogeneous
  # Dirichlet boundary conditions (BCs).
  w = @. 1 + (exp(1) - 1) * x

  # This is still the same RHS since `w` is harmonic! In general, you would need
  # to modify the RHS.
  f = @. -exp(x)

  # Remember to add `w` after computing the solution with homogeneous BCs
  u_homogeneous = A \ f
  uh = u_homogeneous + w
  return x, uh
end


"""
    exercise01_3(internal_solve, analytical_solution)

Perform the task of exercise 1.3a/b and print numerical errors as well as the
experimental order of convergence (EOC) to the screen.
"""
function exercise01_3(internal_solve, analytical_solution)
  nx_values = 10 .* 2 .^ (0:8)

  # setup arrays to store the results
  errors = zeros(length(nx_values))
  order_of_convergence = similar(errors)

  # compute errors in the maximum norm
  for (idx, nx) in enumerate(nx_values)
    x, uh = internal_solve(nx)
    u_ana = analytical_solution.(x)
    errors[idx] = maximum(abs, uh - u_ana)
  end

  # compute experimental order of convergence
  order_of_convergence[1] = NaN # no EOC defined for the first grid
  for idx in 2:length(errors)
    order_of_convergence[idx] = -( log(errors[idx] / errors[idx - 1]) /
                                   log(nx_values[idx] / nx_values[idx - 1]) )
  end

  # print results using the Julia standard library
  println("  nx |   error  |  EOC")
  println("----------------------")
  for idx in eachindex(errors)
    if idx == firstindex(errors)
      @printf("%4d | %.2e | \n", nx_values[idx], errors[idx])
    else
      @printf("%4d | %.2e | %.2f\n", nx_values[idx], errors[idx], order_of_convergence[idx])
    end
  end

  # now do the same but use the package PrettyTables.jl to also get LaTeX code
  # used for the solution hints
  data = hcat(nx_values, errors, order_of_convergence)
  header = ["nx", "error", "EOC"]
  kwargs = (; header, formatters=(ft_printf("%d", 1), ft_printf("%.2e", 2), ft_printf("%.2f", 3)))
  pretty_table(data; kwargs...)
  pretty_table(data; kwargs..., backend=Val(:latex))

  return nothing
end

# If you want to visualize the results, you can run something like
#   using Plots
#   x, uh = solve_exercise01_3a(50); plot(x, uh)


exercise01_3a() = exercise01_3(solve_exercise01_3a, sinpi)
exercise01_3b() = exercise01_3(solve_exercise01_3b, exp)

