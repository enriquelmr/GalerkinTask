using Gridap
using GridapGmsh
using Gridap.Io

function Gridap.DiscreteModelFromFile(filename::AbstractString,::Val{:msh})
  model = GmshDiscreteModel(filename)
  model
end

model = DiscreteModelFromFile("disc_with_hole_3.msh")

order = 1
reffe = ReferenceFE(lagrangian,Float64,order)


V0 = TestFESpace(model,reffe;conformity=:H1,dirichlet_tags=["inner_boundary","outer_boundary"])
Ug = TrialFESpace(V0, [1.0, 2.0]) #the 2 boundary conditions

degree = 2
f(x) = 1.0 #
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)
a(u,v) = ∫( ∇(v)⋅∇(u) )dΩ
b(v) = ∫( v*f )*dΩ


op = AffineFEOperator(a,l,Ug,V0)
u_h = solve(op)
writevtk(Ω,"laplace_task3",cellfields=["u_h"=>u_h])




