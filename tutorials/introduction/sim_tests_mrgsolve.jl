using PuMaS, Plots, LinearAlgebra
using Random
Random.seed!(1234)
choose_covariates() = (isPM = rand(["yes", "no"]),
                    Wt = rand(55:80))

function generate_population(events,nsubs=24)
  pop = Population(map(i -> Subject(id=i,evs=events,cvs=choose_covariates()),1:nsubs))
  return pop
end

m_diffeq = @model begin
    @param   begin
        θ ∈ VectorDomain(7, lower=zeros(7), init=ones(7))
        Ω ∈ PDiagDomain(2)
        σ_prop ∈ RealDomain(init=0.1)
    end

    @random begin
      η ~ MvNormal(Ω)
    end

    @pre begin
        Ka = θ[1]
        TVCL = isPM == "yes" ? θ[2] : θ[2]*θ[6]
        CL = TVCL*exp(η[1])
        V  = θ[3]*exp(η[2])
        lags = [0,θ[4]]
        bioav = [1,θ[5]]
        duration = [0,θ[7]]
    end

    @covariates isPM Wt

     #@dynamics ImmediateAbsorptionModel
      @dynamics begin
          Depot'   = -Ka*Depot
          Central' =  Ka*Depot - (CL/V)*Central
      end

    @derived begin
        cp = @. 1000*(Central / V)
        dv ~ @. Normal(cp, sqrt(cp^2*σ_prop))
    end
end


# Bolus with additional
ev = DosageRegimen(100, ii=24, addl=3)
ev1 =  generate_population(ev)

p = (  θ = [1.5,  #Ka
           1.1,  #CL
           20.0,  #V
           eps(), # lags2
           1, #Bioav
           0.5, # isPM CL
           eps() # duration
           ],
      Ω = Diagonal([0.04,0.04]),
      σ_prop = 0.04
  )


sim1 = simobs(m_diffeq, ev1, p; abstol=1e-14, reltol=1e-14)
simdf1 = DataFrame(sim1)
data1 = read_pumas(simdf1,cvs=[:isPM,:Wt])
res1 = fit(m_diffeq,data1,p,PuMaS.FOCEI())
plot(sim1)
# vpcdf = vpc(m_diffeq,ev1,p,200)
# vpcdf |> plot

# Bolus with lag time and bioav lag2=12.13, bioav=2.23
# FIXME BUG IN lags when more than one event #364
ev = DosageRegimen(100, ii = 24, addl = 3, cmt = 2)
ev2 =  generate_population(ev)

p = (  θ = [1.5,  #Ka
           1.1,  #CL
           20.0,  #V
           12.13, # lags2
           2.23, #Bioav
           0.5, # isPM CL
           0 # duration
           ],
      Ω = Diagonal([0.04,0.04]),
      σ_prop = 0.04
  )


sim2 = simobs(m_diffeq, ev2, p; abstol=1e-14, reltol=1e-14)
plot(sim2)
simdf2 = DataFrame(sim2)
data2 = read_pumas(simdf2,cvs=[:isPM,:Wt])
res2 = fit(m_diffeq,data2,p,PuMaS.FOCEI())

# Infusion with additional
ev = DosageRegimen(100, ii = 24, addl = 3, rate = 100/10, cmt = 2)
ev3 =  generate_population(ev)

p = (  θ = [1.5,  #Ka
           1.1,  #CL
           20.0,  #V
           eps(), # lags2
           1, #Bioav
           0.5, # isPM CL
           eps() # duration
           ],
      Ω = Diagonal([0.04,0.04]),
      σ_prop = 0.04
  )


sim3 = simobs(m_diffeq, ev3, p; abstol=1e-14, reltol=1e-14)
simdf3 = DataFrame(sim3)
data3 = read_pumas(simdf3,cvs=[:isPM,:Wt])
res3 = fit(m_diffeq,data3,p,PuMaS.FOCEI())
plot(sim)


# Infusion doses to depot, with additional
ev = DosageRegimen(100, ii = 24, addl = 3, rate = 100/12, cmt = 1)
ev4 =  generate_population(ev)

p = (  θ = [1.5,  #Ka
           1.1,  #CL
           20.0,  #V
           eps(), # lags2
           1, #Bioav
           eps(), # isPM CL
           0 # duration
           ],
      Ω = PDMat(diagm(0 => [0.04,0.04])),
      σ_prop = 0.00
  )


sim = simobs(m_diffeq, ev4, p; abstol=1e-14, reltol=1e-14)
plot(sim)

# Infusion doses, with additional and lag time = 4.15
# FIXME BUG IN lags when more than one event #364
ev = DosageRegimen(100, ii = 24, addl=3, rate = 100/10, cmt = 2)
ev5 =  generate_population(ev)

p = (  θ = [1.5,  #Ka
           1.1,  #CL
           20.0,  #V
           4.15, # lags2
           1, #Bioav
           0.5, # isPM CL
           0 # duration
           ],
      Ω = PDMat(diagm(0 => [0.04,0.04])),
      σ_prop = 0.00
  )


sim = simobs(m_diffeq, ev5, p; abstol=1e-14, reltol=1e-14)
plot(sim)

# Infusion doses, with lag time = 3.15 and bioav factor =0.412
# FIXME BUG IN lags when more than one event #364
ev = DosageRegimen(100, ii = 24, addl = 3, rate = 100/10, cmt = 2)
ev6 =  generate_population(ev)


p = (  θ = [1.5,  #Ka
           1.1,  #CL
           20.0,  #V
           3.15, # lags2
           0.412, #Bioav
           0.5, # isPM CL
           0 # duration
           ],
      Ω = PDMat(diagm(0 => [0.04,0.04])),
      σ_prop = 0.00
  )


sim = simobs(m_diffeq, ev6, p; abstol=1e-14, reltol=1e-14)
plot(sim)

# Infusion doses, with lag time and bioav factor
ev = DosageRegimen(100, ii = 24, addl = 3, rate = 100/10, ss = 1, cmt = 2)
ev7 =  generate_population(ev)


p = (  θ = [1.5,  #Ka
           1.1,  #CL
           20.0,  #V
           3.16, # lags2
           0.412, #Bioav
           0.5, # isPM CL
           0 # duration
           ],
      Ω = PDMat(diagm(0 => [0.04,0.04])),
      σ_prop = 0.00
  )


sim = simobs(m_diffeq, ev7, p; abstol=1e-14, reltol=1e-14)
plot(sim)

# Infusion doses at steady-state, with bioav factor=0.812
ev = DosageRegimen(100, ii = 12, addl = 4, rate = 100/50, ss = 1, cmt = 2)
ev8 =  generate_population(ev)

p = (  θ = [1.5,  #Ka
           1.1,  #CL
           20.0,  #V
           0, # lags2
           0.812, #Bioav
           0.5, # isPM CL
           0 # duration
           ],
      Ω = PDMat(diagm(0 => [0.04,0.04])),
      σ_prop = 0.00
  )


sim = simobs(m_diffeq, ev8, p; abstol=1e-14, reltol=1e-14)
plot(sim)

# Infusion doses, at steady state
ev = DosageRegimen(100, ii = 12, addl = 3, rate = 100/50, ss = 1, cmt = 2)
ev9 =  generate_population(ev)

p = (  θ = [1.5,  #Ka
           1.1,  #CL
           20.0,  #V
           0, # lags2
           1, #Bioav
           0.5, # isPM CL
           0 # duration
           ],
      Ω = PDMat(diagm(0 => [0.04,0.04])),
      σ_prop = 0.00
  )


sim = simobs(m_diffeq, ev9, p; abstol=1e-14, reltol=1e-14)
plot(sim)

# Infusion doses at steady state, II < DUR, no bioav factor
ev = DosageRegimen(100, ii = 6, addl = 12, rate = round(100/12,digits=5), ss = 1, cmt = 2)
ev10 =  generate_population(ev)

p = (  θ = [1.5,  #Ka
           1.1,  #CL
           20.0,  #V
           0, # lags2
           1, #Bioav
           0.5, # isPM CL
           0 # duration
           ],
      Ω = PDMat(diagm(0 => [0.04,0.04])),
      σ_prop = 0.00
  )


sim = simobs(m_diffeq, ev10, p; abstol=1e-14, reltol=1e-14)
plot(sim)

# Infusion doses at steady state where II == DUR, with bioav factor
ev = DosageRegimen(100, ii = 10, addl = 8, rate = 0.412*100/10,  ss = 1, cmt = 2)
ev11 =  generate_population(ev)

p = (  θ = [1.5,  #Ka
           1.1,  #CL
           20.0,  #V
           0, # lags2
           0.412, #Bioav
           0.5, # isPM CL
           0 # duration
           ],
      Ω = PDMat(diagm(0 => [0.04,0.04])),
      σ_prop = 0.00
  )


sim = simobs(m_diffeq, ev11, p; abstol=1e-14, reltol=1e-14)
plot(sim)


# Infusion doses at steady state, where II == DUR
ev = DosageRegimen(100, ii = 10, addl = 8, rate = 100/10, ss = 1, cmt = 2)
ev12 =  generate_population(ev)

p = (  θ = [1.5,  #Ka
           1.1,  #CL
           20.0,  #V
           0, # lags2
           1, #Bioav
           0.5, # isPM CL
           0 # duration
           ],
      Ω = PDMat(diagm(0 => [0.04,0.04])),
      σ_prop = 0.00
  )


sim = simobs(m_diffeq, ev12, p; abstol=1e-14, reltol=1e-14)
plot(sim)


# Bolus doses at steady state, with bioav factor and lag time
#FIXME lags don't work on subsequent doses
ev = DosageRegimen(100, ii = 24, addl = 3,  ss = 1, cmt = 2)
ev13 =  generate_population(ev)
p = (  θ = [1.5,  #Ka
           1.1,  #CL
           20.0,  #V
           4, # lags2
           0.412, #Bioav
           0.5, # isPM CL
           0 # duration
           ],
      Ω = PDMat(diagm(0 => [0.04,0.04])),
      σ_prop = 0.00
  )


sim = simobs(m_diffeq, ev13, p; abstol=1e-14, reltol=1e-14)
plot(sim)

# Bolus doses with lag time and bioavability factor
ev = DosageRegimen(100, ii = 24, addl = 3, cmt = 2)
ev14 =  generate_population(ev)

p = (  θ = [1.5,  #Ka
           1.1,  #CL
           20.0,  #V
           5, # lags2
           0.412, #Bioav
           0.5, # isPM CL
           0 # duration
           ],
      Ω = PDMat(diagm(0 => [0.04,0.04])),
      σ_prop = 0.00
  )


sim = simobs(m_diffeq, ev14, p; abstol=1e-14, reltol=1e-14)
plot(sim)

# Bolus then infusion
bolus = DosageRegimen(100, cmt = 2,time = 0)
infusion = DosageRegimen(50, time = 13, ii = 24, addl = 2, rate = 24)
ev = DosageRegimen(bolus,infusion)
ev15 =  generate_population(ev)

p = (  θ = [1.5,  #Ka
           1.1,  #CL
           20.0,  #V
           1, # lags2
           1, #Bioav
           0.5, # isPM CL
           0 # duration
           ],
      Ω = PDMat(diagm(0 => [0.04,0.04])),
      σ_prop = 0.00
  )


sim = simobs(m_diffeq, ev15, p; abstol=1e-14, reltol=1e-14)
plot(sim)

# Infusion with modeled duration=9, lag time=5, and bioav factor=0.61
#FIXME DosageRegimen should allow rate value of -2 to allow modeling durations
# rate = 0 allows modeling durations in pumas. Rate=-2 is nonmem convention
ev = DosageRegimen(100, rate = -2, cmt = 2, addl=3,ii=12)
ev16 =  generate_population(ev)

p = (  θ = [1.5,  #Ka
           1.1,  #CL
           20.0,  #V
           5, # lags2
           0.61, #Bioav
           0.5, # isPM CL
           9 # duration
           ],
      Ω = PDMat(diagm(0 => [0.04,0.04])),
      σ_prop = 0.00
  )


sim = simobs(m_diffeq, ev16, p; abstol=1e-14, reltol=1e-14)
plot(sim)

# Infusion with modeled duration=9, at steady state with bioav factor=0.61
#FIXME DosageRegimen should allow rate value of -2 to allow modeling durations
ev = DosageRegimen(100, rate = 0, cmt = 2, ii = 24, addl = 3, ss = 1)
ev17 =  generate_population(ev)


p = (  θ = [1.5,  #Ka
           1.1,  #CL
           20.0,  #V
           0, # lags2
           0.61, #Bioav
           0.5, # isPM CL
           9 # duration
           ],
      Ω = PDMat(diagm(0 => [0.04,0.04])),
      σ_prop = 0.00
  )


sim = simobs(m_diffeq, ev17, p; abstol=1e-14, reltol=1e-14)
plot(sim)


# Reset and dose (EVID 4) with additional
#FIXME does not work as intended
firstdose = DosageRegimen(100, ii = 12, addl = 2, rate = 50)
#FIXME second dose BIOAV should be 0.5 how do we specify this?
seconddose = DosageRegimen(120, evid = 4, time = 50, ii = 12, addl = 3)
ev = DosageRegimen(firstdose,seconddose)
ev18 =  generate_population(ev)

m_diffeq = @model begin
    @param   begin
        θ ∈ VectorDomain(8, lower=zeros(8), init=ones(8))
        Ω ∈ PSDDomain(2)
        σ_prop ∈ RealDomain(init=0.1)
    end

    @random begin
      η ~ MvNormal(Ω)
    end

    @pre begin
        Ka = θ[1]
        TVCL = isPM == "yes" ? θ[2] : θ[2]*θ[6]
        CL = TVCL*exp(η[1])
        V  = θ[3]*exp(η[2])
        lags = [0,θ[4]]
        bioav = time >= 50 ? [1,θ[5]] : [1,θ[8]]
        #bioav = [1,θ[5]]
        duration = [0,θ[7]]
    end

    @covariates isPM Wt

    @dynamics begin
        Depot'   = -Ka*Depot
        Central' =  Ka*Depot - (CL/V)*Central
    end

    @derived begin
        cp = @. 1000*(Central / V)
    end
end


p = (  θ = [1.5,  #Ka
           1.1,  #CL
           20.0,  #V
           0, # lags2
           0.61, #Bioav central
           0.5, # isPM CL
           9, # duration
           0.5, #Bioav central after time=50
           ],
      Ω = PDMat(diagm(0 => [0.04,0.04])),
      σ_prop = 0.00
  )


sim = simobs(m_diffeq, ev18, p; abstol=1e-14, reltol=1e-14)
plot(sim)

# Reset and dose (EVID 4) with additional
#FIXME evid=3 error #371 fixed however assertion rules do not allow creation of population
firstdose = DosageRegimen(100, ii = 12, addl = 3, rate = 50) #bioav here is , BIOAV = 0.61
seconddose = DosageRegimen(0, time = 50, evid = 3, cmt = 2)
thirddose = DosageRegimen(120, time = 54, ii = 16, addl = 2)
#ev = DosageRegimen(firstdose,seconddose,thirddose)
ev = DosageRegimen(DosageRegimen(firstdose,seconddose),thirddose)
ev19 =  generate_population(ev)

p = (  θ = [1.5,  #Ka
           1.1,  #CL
           20.0,  #V
           0, # lags2
           1, #Bioav central
           0.5, # isPM CL
           9, # duration
           0.5, #Bioav central after time=50
           ],
      Ω = PDMat(diagm(0 => [0.04,0.04])),
      σ_prop = 0.00
  )


sim = simobs(m_diffeq, ev19, p; abstol=1e-14, reltol=1e-14)
plot(sim)

# Steady state 1 and 2
firstdose = DosageRegimen(100, ii = 24, addl = 3, ss = 1)
seconddose = DosageRegimen(50, time = 12, ii = 24, addl = 3, ss = 2)
ev = DosageRegimen(firstdose,seconddose)
ev20 =  generate_population(ev)

p = (  θ = [1.5,  #Ka
           1.1,  #CL
           20.0,  #V
           0, # lags2
           1, #Bioav
           0.5, # isPM CL
           0 # duration
           ],
      Ω = PDMat(diagm(0 => [0.04,0.04])),
      σ_prop = 0.00
  )


sim = simobs(m_diffeq, ev20, p; abstol=1e-14, reltol=1e-14)
plot(sim)