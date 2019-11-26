using Pumas, Test

@test Central1().syms == (:Central,)
@test Pumas.DiffEqBase.has_syms(Central1())
@test Depots1Central1().syms == (:Depot, :Central)
@test Pumas.DiffEqBase.has_syms(Depots1Central1())
@test Depots2Central1().syms == (:Depot1, :Depot2, :Central)
@test Pumas.DiffEqBase.has_syms(Depots2Central1())

#==
  Central1

  Eigen values and vectors are very simple: they're just -CL/V and [1] respectively.
==#
p = (CL=rand(), V=rand())
ocm = Pumas.LinearAlgebra.eigen(Central1(), p)
@test all(ocm[1] .== -p.CL/p.V)
@test all(ocm[2] .== [1.0])

#==
  Depots1Central1

  Eigen values are just: Λ = [-Ka, -CL/V]
  Eigen vectors are [Λ[1]/Λ[2] - 1, 0] and [0,1]
==#

p = (Ka=rand(), CL=rand(), V=rand())
ocm = Pumas.LinearAlgebra.eigen(Depots1Central1(), p)
@test all(ocm[1] .== [-p.Ka, -p.CL/p.V])
@test all(ocm[2] .== [ocm[1][2]/ocm[1][1]-1  0.0; 1.0 1.0])

#==
  Two compartment, one Central one Peripheral. Compare against solution of nu-
  merical matrix given to eigen solver.
==#

p = (CL=0.1, Vc=1.0, Vp=2.0, Q=0.5)
ocm = Pumas.LinearAlgebra.eigen(Central1Periph1(), p)
λ1 = -(17+sqrt(249))/40
λ2 = -(17-sqrt(249))/40
@test all(ocm[1] .≈ [λ1, λ2])
@test all(ocm[2] .≈ [2*λ1+1/2 2*λ2+1/2; 1 1])
ocm[2]
p = (CL=0.1, Vc=5.0, Vp=2.0, Q=0.5)
ocm = Pumas.LinearAlgebra.eigen(Central1Periph1(), p)
λ1 = -(37+sqrt(1169))/200
λ2 = -(37-sqrt(1169))/200
v1 = -(-13+sqrt(1169))/20
v2 = -(-13-sqrt(1169))/20
@test all(ocm[1] .≈ [λ1, λ2])
@test all(ocm[2] .≈ [v1 v2; 1 1])

#==
  Depots1Central1Periph1
==#

p = (Ka=0.05, CL=0.1, Vc=1.0, Vp=2.0, Q=0.5)
ocm = Pumas.LinearAlgebra.eigen(Depots1Central1Periph1(), p)
λ1 = -(17+sqrt(249))/40
λ2 = -(17-sqrt(249))/40
@test all(ocm[1] .≈ [-p.Ka, λ1, λ2])
@test all(ocm[2] .≈ [-3/5 0 0; 2/5 2*λ1+1/2 2*λ2+1/2; 1 1 1])

p = (Ka=0.05, CL=0.1, Vc=5.0, Vp=2.0, Q=0.5)
ocm = Pumas.LinearAlgebra.eigen(Depots1Central1Periph1(), p)
λ1 = -(37+sqrt(1169))/200
λ2 = -(37-sqrt(1169))/200
v1 = -(-13+sqrt(1169))/20
v2 = -(-13-sqrt(1169))/20
@test all(ocm[1] .≈ [-p.Ka, λ1, λ2])
@test all(ocm[2] .≈ [-2.2 0 0; 2.0 v1 v2; 1 1 1])

#==
  Central1Periph1MetaPeriph1
==#
p = (CL1=0.01, CL2=0.1, V1=5.0, Vp1=2.0, V2=6.0, Vp2=3.0, Q1=0.5, Q2=1.2, T=1.0)
ocm = Pumas.LinearAlgebra.eigen(Central1Periph1MetaPeriph1(), p)
λ1 = -(138+7*sqrt(131))/500
λ2 = -(138-7*sqrt(131))/500
λ3 = -(37+sqrt(1273))/120
λ4 = -(37-sqrt(1273))/120
V = [-(46261+679*sqrt(131))/30000 -(46261-679*sqrt(131))/30000 0 0;(17+252*sqrt(131))/3000 (17-252*sqrt(131))/3000  0  0; (62-7*sqrt(131))/100    (62+7*sqrt(131))/100     (11-sqrt(1273))/24 (11+sqrt(1273))/24;1 1 1 1]
@test all(ocm[1] .≈ [λ1, λ2, λ3, λ4])
@test all(ocm[2] .≈ V)
ocm[2]
ocm[1]
