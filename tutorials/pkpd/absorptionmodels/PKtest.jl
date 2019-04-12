
model = @model begin
  @param begin
    θ ∈ VectorDomain(3)
    Ω ∈ PSDDomain(3)
    σ_prop ∈ RealDomain(init=0.1)
  end

  @random begin
    η ~ MvNormal(Ω)
  end

  @covariates Wt amt

  @pre begin
    CL = θ[1]*(Wt/70)^0.75*exp(η[1])
    V  = θ[2]*(Wt/70)*exp(η[2])
    Ktr = θ[3]*exp(η[3])
    K  = CL/V
    N = 5 # number of compartments
  end

  @dynamics begin
    Central' = Ktr*amt*((Ktr*t)^N/factorial(N))*exp(-Ktr*t)-Central*K
  end

  @vars begin
    conc = Central/V
  end

  @derived begin
    dv ~ @.Normal(conc,sqrt(conc^2*σ_prop + eps()))
  end

end


## transit


model = @model begin
  @param begin
    θ ∈ VectorDomain(6)
    Ω ∈ PSDDomain(6)
    σ_prop ∈ RealDomain(init=0.1)
  end

  @random begin
    η ~ MvNormal(Ω)
  end

  @covariates Wt amt

  @pre begin
    CL = θ[1]*(Wt/70)^0.75*exp(η[1])
    V  = θ[2]*(Wt/70)*exp(η[2])
    BIO = θ[3]*exp(η[3]) # bioavailability
    Ka = θ[4]*exp(η[4])
    MTT =θ[5]*exp(η[5])
    NN = θ[6]*exp(η[6])
    K  = CL/V
    Ktr = (NN+1)/MTT
    LNFAC=log(2.5066)+(NN+0.5)*log(NN)-NN
  end

  @dynamics begin
    Transit' = exp(log(BIO*amt+.00001)+log(Ktr)+NN*log(Ktr*t+0.00001)-Ktr*t-LNFAC)-Ka*Transit
    Central' = Ka*Transit-Central*K
  end

  @vars begin
    conc = Central/V
  end

  @derived begin
    dv ~ @.Normal(conc,sqrt(conc^2*σ_prop + eps()))
  end

end


model = @model begin

  @param begin
    θ ∈ VectorDomain(5)
    Ω ∈ PSDDomain(3)
    σ_prop ∈ RealDomain(init=0.1)
  end

  @random begin
    η ~ MvNormal(Ω)
  end

  @covariates Wt amt

  @pre begin
    CL = θ[1]*(Wt/70)^0.75*exp(η[1])
    V  = θ[2]*(Wt/70)*exp(η[2])
    Ka = θ[3]
    duration = [0, θ[4]]
    bioav = [1-θ[5], θ[5]]
  end

  @dynamics begin
    Depot'   = -Ka*Depot
    Central' =  Ka*Depot - Central*CL/V
  end

  @vars begin
    conc = Central/V
  end

  @derived begin
    dv ~ @.Normal(conc,sqrt(conc^2*σ_prop + eps()))
  end

end
