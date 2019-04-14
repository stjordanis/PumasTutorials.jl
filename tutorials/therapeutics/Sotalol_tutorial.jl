---
title : Sotalol dosing simuations: Oral and IV to Oral conversion
author : Elyes Dahmane, Vijay Ivaturi
date :  April 14, 2019
---

```julia
using PuMaS
using LinearAlgebra
using Plots
```

```julia
#using Random
#How to set a seed ? Random.seed!(04152019)
```

# Introduction

The objectives of this tutorial is to demonstrate with an example how to perform PKPD simulations of different dosing regimen and how to plots the PK/PD profile in Julia,
using the package PuMaS, and its dependencies eg. the package LinearAlgebra (what specific functions ?) and Plots.

The 1st step will be the creation of different sotalol dosing regimen to simulate.
The 2nd step will be the simulation of the mean profile for each dosing regimen an plotting the PK/PD profiles.
The 3rd step will be to simulate a population.

# Building the dataset for the different Oral and IV to Oral vonversion dosing regimen for sotalol

Creation of 12 different dosing regimen:

Sotalol is usually orally administred at 3 dose levels: 80 mg, 120 mg or 160 mg twice daily (BID) in patient with normal to mild renal function (> 60 mL/min).
In patient with moderate renal impairement (40 - 60 mL/min), sotalol is administered once daily QD.

We are going to simulate the sotalol concentrations under each dosage:
reg_PO_80, reg _PO_120 and reg _PO_160 are regimen were the drug administred PO every 12 hours for three days (6 doses in total) under a dose of eith 80 mg, 120 mg or 160 mg.

reg_PO_80_CKI, reg _PO_120_CKI and reg _PO_160_CKI are regimen were the drug administred PO every 24 hours for three days (6 doses in total) under a dose of eith 80 mg, 120 mg or 160 mg.

```julia
##caveat: inconvenient that I cannot specify the name of first argument ... at least amt=[40,80,80]
## 80 mg BID and QD dosing regimen:
reg_IVPO_80 = DosageRegimen([40,80,80], cmt=[2,1,1], time=[0,1,12],
                        ii=[0,0,12], addl=[0,0,4], rate=[40/1,0,0])
reg_IVPO_80_CKI = DosageRegimen([40,80,80], cmt=[2,1,1], time=[0,1,24],
                        ii=[0,0,24], addl=[0,0,4], rate=[40/1,0,0])

reg_PO_80 = DosageRegimen(80, cmt=1, time=0, ii=12, addl=5, rate=0)
reg_PO_80_CKI = DosageRegimen(80, cmt=1, time=0, ii=24, addl=5, rate=0)

## 120 mg BID and QD dosing regimen:
reg_IVPO_120 = DosageRegimen([40,20,120,120], cmt=[2,2,1,1], time=[0,1,1.5,12],
                        ii=[0,0,0,12], addl=[0,0,0,4], rate=[40/1,20/0.5,0,0])
reg_IVPO_120_CKI = DosageRegimen([40,10,120,120], cmt=[2,2,1,1], time=[0,1,1.5,24],
                        ii=[0,0,0,24], addl=[0,0,0,4], rate=[40/1,10/0.5,0,0])

reg_PO_120 = DosageRegimen(120, cmt=1, time=0, ii=12, addl=5, rate=0)
reg_PO_120_CKI = DosageRegimen(120, cmt=1, time=0, ii=24, addl=5, rate=0)

## 160 mg BID and QD dosing regimen:
reg_IVPO_160 = DosageRegimen([40,20,20,160,160], cmt=[2,2,2,1,1], time=[0,1,1.5,2,12],
                        ii=[0,0,0,0,12], addl=[0,0,0,0,4], rate=[40/1,20/0.5,20/0.5,0,0])
reg_IVPO_160_CKI = DosageRegimen([40,10,20,160,160], cmt=[2,2,2,1,1], time=[0,1,1.5,2,24],
                        ii=[0,0,0,0,24], addl=[0,0,0,0,4], rate=[40/1,10/0.5,20/0.5,0,0])

reg_PO_160 = DosageRegimen(160, cmt=1, time=0, ii=12, addl=5, rate=0)
reg_PO_160_CKI = DosageRegimen(160, cmt=1, time=0, ii=24, addl=5, rate=0)
```

#Creation of 12 typical subjects (1 subject per dosing regimen)
##evs and cvs are argument names for events and covariates.
s1 = Subject(id=1, evs=reg_IVPO_80, cvs=(BW = 70, CKI=0, dose=["iv1","po","po"], reg="IV_PO_80mg"))
s2 = Subject(id=2, evs=reg_PO_80, cvs = (BW = 70, CKI=0))
s3 = Subject(id=3, evs=reg_IVPO_120, cvs=(BW = 70, CKI=0, dose=["iv1","iv2","po","po"]))
s4 = Subject(id=4, evs=reg_PO_120, cvs=(BW = 70, CKI=0))
s5 = Subject(id=5, evs=reg_IVPO_160, cvs=(BW = 70, CKI=0, dose=["iv1","iv2","iv3","po","po"]))
s6 = Subject(id=6, evs=reg_PO_160, cvs=(BW = 70, CKI=0))

s7 = Subject(id=7, evs=reg_IVPO_80_CKI, cvs=(BW = 70, CKI=1, dose=["iv1","po","po"]))
s8 = Subject(id=8, evs=reg_PO_80_CKI, cvs=(BW = 70,CKI=1))
s9 = Subject(id=9, evs=reg_IVPO_120_CKI, cvs=(BW = 70, CKI=1, dose=["iv1","iv2","po","po"]))
s10 = Subject(id=10, evs=reg_PO_120_CKI, cvs=(BW = 70, CKI=1))
s11 = Subject(id=11, evs=reg_IVPO_160_CKI, cvs=(BW = 70, CKI=1, dose=["iv1","iv2","iv3","po","po"]))
s12 = Subject(id=12, evs=reg_PO_160_CKI, cvs=(BW = 70, CKI=1))

#Creation of one big dataset by merging the 12 typical subjects (1 subject per dosing regimen)
mypop=Population([s1,s2,s3,s4,s5,s6, s7,s8,s9,s10,s11,s12])

#mypop1=Population([s1])
#mypop2=Population([s2])

#Creation of the PKPD model for simulation of the mean PK/PD profiles for each dosing regimen
## what is the meaning of @pre; Answer: abreviation of pre-processing.
## is bioav and lags proprietary names
## what is the meaning of the first value in the array of bioav and lags.
## what is the meaning of apostrophe in @dynamics
## what should I do if I want to specify K20 or Kel=CL/Vc. shoud it be done in @vars ?

##PK model

sotalol_PKmodel = @model begin

    @param   begin
        θ ∈ VectorDomain(7, lower=zeros(7), init=ones(7))
        σ_prop ∈ RealDomain(init=0.1)
    end

    @random begin
      η ~ MvNormal(Matrix{Float64}(I, 7, 7))
    end

    @pre begin
        bioav = [1, θ[1]]#*exp(η[1])]
        Ka = θ[2]#*exp(η[2])
        lags = [θ[3],0]#*exp(η[3])]
        CL = θ[4]*(BW/70)^0.75*(0.5)^CKI#*exp(η[4])
        Vc  = θ[5]*(BW/70)#*exp(η[5])
        Q  = θ[6]*(BW/70)^0.75#*exp(η[6])
        Vp  = θ[7]*(BW/70)#*exp(η[7])
    end

    @covariates BW CKI

    @dynamics begin
        Depot'   = -Ka*Depot
        Central' =  Ka*Depot - (CL/Vc)*Central - (Q/Vc)*Central + (Q/Vp)*Peripheral
        Peripheral' = (Q/Vc)*Central - (Q/Vp)*Peripheral
    end

    @derived begin
        CP = @. (Central / (Vc/1000))
        DV = @. Normal(CP, sqrt(CP^2*σ_prop))
    end

end


parameters_PK = ( θ = [1.0667, #bioav
                        0.605, #Ka
                        0.231, #alag1
                        12, #CL
                        77.1, #Vc
                        9.22, #Q
                        52.3 #Vp
                        ],
                    #Ω = [0,0,0,0,0,0,0],
                    σ_prop = 0)


sim2=simobs(sotalol_PKmodel, s2, parameters_PK)
plot(sim2)

#Plots shows the title population simulation eventhough it is an individual simulation.
#cannot control or visualize what is outputed in sim. (total black box)

##PK/PD model

## Why do I need +eps( ) in the DV calculation
## How to specify 2 different residual errors
## Because σ_prop is not specified as variance from a matrix, do I need to specify a lower boundary
## can we define or update BW in parameters_PK ?

sotalol_PKPD_model = @model begin

    @param   begin
        θ ∈ VectorDomain(9, lower=zeros(9), init=ones(9))
        Ω ∈ PSDDomain(9)
        σ_prop_PK ∈ RealDomain(init=0.1)
        σ_prop_PD ∈ RealDomain(init=0.1)
    end

    @random begin
      η ~ MvNormal(Ω)
    end

    @pre begin
        bioav = [1,θ[1]*exp(η[1])]
        Ka = θ[2]*exp(η[2])
        lags = [θ[3]*exp(η[3]),0]
        CL = θ[4]*(BW/70)^0.75*(0.5)^CKI*exp(η[4])
        Vc  = θ[5]*(BW/70)*exp(η[5])
        Q  = θ[6]*(BW/70)^0.75*exp(η[6])
        Vp  = θ[7]*(BW/70)*exp(η[7])

        E0 = θ[8]*exp(η[8])
        Slope = θ[9]*exp(η[9])
    end

    @covariates BW CKI

    @dynamics begin
        Depot'   = -Ka*Depot
        Central' =  Ka*Depot - (CL/Vc)*Central - (Q/Vc)*Central + (Q/Vp)*Peripheral
        Peripheral' = (Q/Vc)*Central - (Q/Vp)*Peripheral
    end

    @derived begin
        CP = @. (Central / (Vc/1000))
        DV = @. Normal(CP, sqrt(CP^2*σ_prop_PK))

        E = @. E0 + (Slope*CP)
        QT = @. Normal(E, sqrt(σ_prop_PD))
    end

end


parameters_PKPD = ( θ = [1.0667, #bioav
                            0.605, #Ka
                            0.231, #alag1
                            12, #CL
                            77.1, #Vc
                            9.22, #Q
                            52.3, #Vp
                            405, #E0
                            0.0158 #Slope
                            ],
                    Ω = PDMat(diagm(0 => [1E-16, 0.4692, 1E-16, 0.0198, 0.0713, 1E-16, 1E-16, 0.0025, 0.034])),
                    #Ω = PDMat(diagm(0 => [1E-16, 1E-16, 1E-16, 1E-16, 1E-16, 1E-16, 1E-16, 1E-16, 1E-16])),
                    σ_prop_PK = 0.06, # made it smaller because it was ridiculously high 0.2285 (%CV of 47.8),
                    σ_prop_PD = 19.1
                    )


sim1=simobs(sotalol_PKPD_model, s1, parameters_PKPD)
plot(sim1)

sim=simobs(sotalol_PKPD_model, mypop, parameters_PKPD)
#sim=simobs(sotalol_PKPD_model, mypop, parameters_PKPD, obstimes=0:24)
plot(sim)

##randeffs = (η = Diagonal([0, 0, 0, 0, 0, 0, 0]),)

#Creation of a population of 1000 subject per dosing regimen
##evs and cvs are argument names for events and covariates.
# to round you have to specify the argument digit=s if the number of digits is different from 0,  round(8.315, digits=1)

pop1 = Population(map(i -> Subject(id=i, evs=reg_IVPO_80, cvs=(BW = round.(rand(Normal(70, 10))), CKI=0, dose=["iv1","po","po"])), 1:1000))
pop2 = Population(map(i -> Subject(id=i, evs=reg_PO_80, cvs = (BW = round.(rand(Normal(70, 10))), CKI=0)), 1001:2000))
pop3 = Population(map(i -> Subject(id=i, evs=reg_IVPO_120, cvs=(BW = round.(rand(Normal(70, 10))), CKI=0, dose=["iv1","iv2","po","po"])), 2001:3000))
pop4 = Population(map(i -> Subject(id=i, evs=reg_PO_120, cvs=(BW = round.(rand(Normal(70, 10))), CKI=0)), 3001:4000))
pop5 = Population(map(i -> Subject(id=i, evs=reg_IVPO_160, cvs=(BW = round.(rand(Normal(70, 10))), CKI=0, dose=["iv1","iv2","iv3","po","po"])), 4001:5000))
pop6 = Population(map(i -> Subject(id=i, evs=reg_PO_160, cvs=(BW = round.(rand(Normal(70, 10))), CKI=0)), 5001:6000))

pop7 = Population(map(i -> Subject(id=i, evs=reg_IVPO_80_CKI, cvs=(BW = round.(rand(Normal(70, 10))), CKI=1, dose=["iv1","po","po"])), 6001:7000))
pop8 = Population(map(i -> Subject(id=i, evs=reg_PO_80_CKI, cvs=(BW = round.(rand(Normal(70, 10))),CKI=1)), 7001:8000))
pop9 = Population(map(i -> Subject(id=i, evs=reg_IVPO_120_CKI, cvs=(BW = round.(rand(Normal(70, 10))), CKI=1, dose=["iv1","iv2","po","po"])), 8001:9000))
pop10 = Population(map(i -> Subject(id=i, evs=reg_PO_120_CKI, cvs=(BW = round.(rand(Normal(70, 10))), CKI=1)), 9001:10000))
pop11 = Population(map(i -> Subject(id=i, evs=reg_IVPO_160_CKI, cvs=(BW = round.(rand(Normal(70, 10))), CKI=1, dose=["iv1","iv2","iv3","po","po"])), 10001:11000))
pop12 = Population(map(i -> Subject(id=i, evs=reg_PO_160_CKI, cvs=(BW =round.(rand(Normal(70, 10))), CKI=1)), 11001:12000))

#Creation of one big dataset by merging the 12 typical subjects (1 subject per dosing regimen)
#merging different populations is different than merging different subjects
#what vcat means ?
mybigpop=Population(vcat(pop1,pop2,pop3,pop4,pop5,pop6, pop7,pop8,pop9,pop10,pop11,pop12))

sim_pop=simobs(sotalol_PKPD_model, mybigpop, parameters_PKPD)

plot(sim_pop)


-
