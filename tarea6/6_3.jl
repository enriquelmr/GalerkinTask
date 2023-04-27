using Pkg
using OrdinaryDiffEq
using Trixi
using Plots
Pkg.activate(@__DIR__)
Pkg.instantiate()



function initial_condition_6(x, t, equations::CompressibleEulerEquations3D)
    
    ρ =1 # density
    v_1=sin(x[1])*cos(x[2])*cos(x[3]) 
    v_2=-cos(x[1])*sin(x[2])*cos(x[3]) 
    v_3=0 
    p=(100/1.4)+((1/16)*(cos(2*x[1])*cos(2*x[3])+2*cos(2*x[2])+2*cos(2*x[1])+cos(2*x[2])*cos(2*x[3]))) #pressure distribution
    
    
    return prim2cons(SVector(ρ,v_1,v_2,v_3,p),equations)
    
end

function ex6_3()
    equations = CompressibleEulerEquations3D(1.4)#define the value of our λ = 1.4, and we got the compresible euler eqautions
    volume_flux = flux_ranocha
    surfacee_flux =flux_lax_friedrichs
    
    solver = DGSEM(polydeg=3, surface_flux=surfacee_flux,
               volume_integral=VolumeIntegralFluxDifferencing(volume_flux))
    coordinates_min=(-π,-π,-π)
    coordinates_max=(π,π,π)
    mesh = TreeMesh(coordinates_min, coordinates_max; periodicity=true, initial_refinement_level=3, n_cells_max=10000)
    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_6, solver)
    
    
    tspan = (0.0, 5.0)
    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()

    analysis_interval = 100
    analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

    alive_callback = AliveCallback(analysis_interval=analysis_interval)

    save_solution = SaveSolutionCallback(interval=100,
                                     save_initial_solution=true,
                                     save_final_solution=true)

    stepsize_callback = StepsizeCallback(cfl=1.0)

    callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        stepsize_callback);
    
    
    sol = solve(ode, RK4(), abstol=1.0e-6, reltol=1.0e-6, save_everystep=false, callback=callbacks);
    
    pl=plot(sol)

    savefig(pl,"sol.png")
    
end
    




