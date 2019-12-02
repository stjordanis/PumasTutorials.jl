export Central1, Depots1Central1, Depots2Central1,
       Central1Periph1, Depots1Central1Periph1 ,
       Central1Periph1Meta1, Central1Periph1MetaPeriph1

abstract type ExplicitModel end

# Generic ExplicitModel solver. Uses an analytical eigen solution.
function _analytical_solve(m::M, t, t₀, amounts, doses, p, rates) where M<:ExplicitModel
  amt₀ = amounts + doses   # initial values for cmt's + new doses
  Λ, 𝕍 = eigen(m, p)

  # We avoid the extra exp calls, but could have written:
  # Dh  = Diagonal(@SVector(exp.(λ * (_t - _t₀)))
  # Dp  = Diagonal(@SVector(expm1.(λ * (_t - _t₀))./λ))
  # We could also have written:
  # Dp = Diagonal(expm1.(Λ * (t - t₀)) ./ Λ)
  # Dh = Dp .* Λ + I
  # but Diagonal{StaticVector} falls back to Array operations. Instead we write:
  dp = expm1.(Λ * (t - t₀)) ./ Λ
  dh = dp .* Λ .+ 1

  # We cannot * here because of Array fallback for Diagonal{StaticVector}
  # amtₜ = 𝕍*(Dp*(𝕍\rates) + Dh*(𝕍\amt₀)) # could derive inverse here
  amtₜ = 𝕍*(dp.*(𝕍\rates) + dh.*(𝕍\amt₀)) # could derive inverse here

  return SLVector(NamedTuple{varnames(M)}(amtₜ))
end
DiffEqBase.has_syms(x::ExplicitModel) = true
Base.getproperty(x::ExplicitModel, symbol::Symbol) = symbol == :syms ? Pumas.varnames(typeof(x)) : getfield(x, symbol)

struct Central1 <: ExplicitModel end
(m::Central1)(args...) = _analytical_solve(m, args...)
@inline function LinearAlgebra.eigen(::Central1, p)
  Ke = p.CL/p.V
  T = typeof(Ke)

  Λ = @SVector([-Ke])
  𝕍 = @SMatrix([T(1)])

  return Λ, 𝕍
end
varnames(::Type{Central1}) = (:Central,)
pk_init(::Central1) = SLVector(Central=0.0)

struct Depots1Central1 <: ExplicitModel end
(m::Depots1Central1)(args...) = _analytical_solve(m, args...)
@inline function LinearAlgebra.eigen(::Depots1Central1, p)
    a = p.Ka
    e = p.CL/p.V

    Λ = @SVector([-a, -e])
    v = e/a - 1
    𝕍 = @SMatrix([v 0;
                  1 1])
    return Λ, 𝕍
end
varnames(::Type{Depots1Central1}) = (:Depot, :Central)
pk_init(::Depots1Central1) = SLVector(Depot=0.0,Central=0.0)

struct Depots2Central1 <: ExplicitModel end
(m::Depots2Central1)(args...) = _analytical_solve(m, args...)
@inline function LinearAlgebra.eigen(::Depots2Central1, p)
    a = p.Ka1
    b = p.Ka2
    e = p.CL/p.V

    frac1 = (e-a)/a
    invfrac1 = inv(frac1)

    frac2 = (e-b)/b
    invfrac2 = inv(frac2)

    Λ = @SVector([-a, -b, -e])

    v1 = -1 + e/a
    v2 = -1 + e/b
    𝕍 = @SMatrix([frac1 0     0;
                  0     frac2 0;
                  1     1     1])

    return Λ, 𝕍
end
varnames(::Type{Depots2Central1}) = (:Depot1, :Depot2, :Central)
pk_init(::Depots2Central1) = SLVector(Depot1=0.0,Depot2=0.0,Central=0.0)

# b is from actual cmt to peri, c is back
struct Central1Periph1 <: ExplicitModel end
_V(::Central1Periph1, Λ, b, c) = @SMatrix([(Λ[1]+c)/b (Λ[2]+c)/b])
function _Λ(::Central1Periph1, a, b, c)
  A = a + b + c
  S = sqrt(A^2-4*a*c)
  Λ = @SVector([-(A+S)/2, -(A-S)/2])
end
(m::Central1Periph1)(args...) = _analytical_solve(m, args...)
@inline function LinearAlgebra.eigen(m::Central1Periph1, p)
    a = p.CL/p.Vc
    b = p.Q/p.Vc
    c = p.Q/p.Vp

    Λ = _Λ(m, a, b, c)
    𝕍 = vcat(_V(m, Λ, b, c), @SMatrix([1 1]))

    return Λ, 𝕍
end
varnames(::Type{Central1Periph1}) = (:Central, :Peripheral)
pk_init(::Central1Periph1) = SLVector(Central=0.0, Peripheral=0.0)

struct Depots1Central1Periph1  <: ExplicitModel end
(m::Depots1Central1Periph1 )(args...) = _analytical_solve(m, args...)
@inline function LinearAlgebra.eigen(::Depots1Central1Periph1 , p)
  k = p.Ka
  a = p.CL/p.Vc
  b = p.Q/p.Vc
  c = p.Q/p.Vp

  A = a + b + c

  Λ, 𝕍 = eigen(Central1Periph1(), p)
  Λ = pushfirst(Λ, -k)

  𝕍 = vcat(@SMatrix([0 0;]), 𝕍) # pad with zeros
  v_depot = @SMatrix([((k-A)+a*c/k)/b; (c-k)/b; 1])
  𝕍 = hcat(v_depot, 𝕍)

  return Λ, 𝕍, inv(𝕍)
end
varnames(::Type{Depots1Central1Periph1 }) = (:Depot, :Central, :Peripheral)
pk_init(::Depots1Central1Periph1 ) = SLVector(Depot=0.0, Central=0.0, Peripheral=0.0)


# use Vc and Vm
struct Central1Periph1MetaPeriph1 <: ExplicitModel end # 011?
(m::Central1Periph1MetaPeriph1)(args...) = _analytical_solve(m, args...)
@inline function LinearAlgebra.eigen(::Central1Periph1MetaPeriph1, p)
  a = p.CL1/p.V1
  b = p.Q1/p.V1
  c = p.Q1/p.Vp1
  d = p.T/p.V1
  e = p.CL2/p.V2
  f = p.Q2/p.V2
  h = p.Q2/p.Vp2

  β = a + b
  ϕ = e + f

  m′ = Central1Periph1()
  Λ = vcat(_Λ(m′, a, b, c),  _Λ(m′, e, f, h))

  v1_3 = ( Λ[1] + h)/f
  v1_1 = ((Λ[1] + ϕ) * v1_3 - h)/d
  v1_2 = ( Λ[1] + β) * (v1_1 + h/d)/c - (Λ[1] + β)*h/(c*d)

  v2_3 = ( Λ[2] + h)/f
  v2_1 = ((Λ[2] + ϕ) * v2_3 - h)/d
  v2_2 = ( Λ[2] + β) * (v2_1 + h/d)/c - (Λ[2] + β)*h/(c*d)


  v3_3 = (Λ[3] + h)/f
  v4_3 = (Λ[4] + h)/f

  𝕍 = @SMatrix([v1_1  v2_1  0   0  ;
                v1_2  v2_2  0   0  ;
                v1_3  v2_3  v3_3 v4_3;
                1     1    1   1])

  return Λ, 𝕍
end
varnames(::Type{Central1Periph1MetaPeriph1}) = (:Central, :CPeripheral, :Metabolite, :MPeripheral)
pk_init(::Central1Periph1MetaPeriph1) = SLVector(Central=0.0, CPeripheral=0.0, Metabolite=0.0, MPeripheral=0.0
)

# use Vc and Vm
struct Central1Periph1Meta1 <: ExplicitModel end # 011?
(m::Central1Periph1Meta1)(args...) = _analytical_solve(m, args...)
@inline function LinearAlgebra.eigen(m::Central1Periph1Meta1, p)
  a = p.CL1/p.V1
  b = p.Q1/p.V1
  c = p.Q1/p.Vp1
  d = p.T/p.V1
  e = p.CL2/p.V2

  β = a + b
  Λ = vcat(_Λ(Central1Periph1(), a, b, c), @SVector([-e]))

  v1_1 = (Λ[1] + e)/d
  v1_2 = (Λ[1] + β)*v1_1/c
  v2_1 = (Λ[2] + e)/d
  v2_2 = (Λ[2] + β)*v2_1/c

  𝕍 = @SMatrix([v1_1 v2_1 0;
                v1_2 v2_2 0;
                1    1    1])

  return Λ, 𝕍
end
varnames(::Type{Central1Periph1Meta1}) = (:Central, :CPeripheral, :Metabolite)
pk_init(::Central1Periph1Meta1) = SLVector(Central=0.0, CPeripheral=0.0, Metabolite=0.0)
