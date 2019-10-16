using Pumas, LinearAlgebra, Plots
using CSV, DataFrames, Query
using GLM, DiffEqSensitivity, Random
using DataStructures, DataFramesMeta, GLM
using Serialization


model = @model begin

    @param begin
        fup ∈ RealDomain(init = 0.42)
        fumic ∈ RealDomain(init = 0.711)
        WEIGHT ∈ RealDomain(init = 73)
        MPPGL ∈ RealDomain(init = 30.3)
        MPPGI ∈ RealDomain(init = 0)
        C_OUTPUT ∈ RealDomain(init = 6.5)
        VmaxH ∈ RealDomain(init = 40)
        VmaxG ∈ RealDomain(init = 40)
        KmH ∈ RealDomain(init = 9.3)
        KmG ∈ RealDomain(init = 9.3)
        BP ∈ RealDomain(init = 1)

        Kpad ∈ RealDomain(init = 9.89)
        Kpbo ∈ RealDomain(init = 7.91)
        Kpbr ∈ RealDomain(init = 7.35)
        Kpgu ∈ RealDomain(init = 5.82)
        Kphe ∈ RealDomain(init = 1.95)
        Kpki ∈ RealDomain(init = 2.9)
        Kpli ∈ RealDomain(init = 4.66)
        Kplu ∈ RealDomain(init = 0.83)
        Kpmu ∈ RealDomain(init = 2.94)
        Kpsp ∈ RealDomain(init = 2.96)
        Kpre ∈ RealDomain(init = 4)

        MW ∈ RealDomain(init = 349.317)
        logP ∈ RealDomain(init = 2.56)
        S_lumen ∈ RealDomain(init = 0.39*1000)
        L ∈ RealDomain(init = 280)
        d ∈ RealDomain(init = 2.5)
        PF ∈ RealDomain(init = 1.57)
        VF ∈ RealDomain(init = 6.5)
        MF ∈ RealDomain(init = 13)
        ITT ∈ RealDomain(init = 3.32)
        A ∈ RealDomain(init = 7440)
        B ∈ RealDomain(init = 1e7)
        alpha ∈ RealDomain(init = 0.6)
        beta ∈ RealDomain(init = 4.395)
        fabs ∈ RealDomain(init = 1)
        fdis ∈ RealDomain(init = 1)
        fperm ∈ RealDomain(init = 1)

        Vad ∈ RealDomain(init = 18.2)
        Vbo ∈ RealDomain(init =10.5)
        Vbr ∈ RealDomain(init =1.45)
        VguWall ∈ RealDomain(init =0.65)
        VguLumen ∈ RealDomain(init =0.35)

        Vhe ∈ RealDomain(init =0.33)
        Vki ∈ RealDomain(init =0.31)
        Vli ∈ RealDomain(init =1.8)
        Vlu ∈ RealDomain(init =0.5)
        Vmu ∈ RealDomain(init =29)
        Vsp ∈ RealDomain(init =0.15)
        Vbl ∈ RealDomain(init =5.6)

        FQad ∈ RealDomain(init = 0.05)
        FQbo ∈ RealDomain(init = 0.05)
        FQbr ∈ RealDomain(init = 0.12)
        FQgu ∈ RealDomain(init = 0.16)
        FQhe ∈ RealDomain(init = 0.04)
        FQki ∈ RealDomain(init = 0.19)
        FQli ∈ RealDomain(init = 0.255)
        FQmu ∈ RealDomain(init = 0.17)
        FQsp ∈ RealDomain(init = 0.03)
    end

    @pre begin

        Vgu = VguWall + VguLumen
        Vve = 0.705*Vbl
        Var = 0.295*Vbl
        Vre = WEIGHT - (Vli+Vki+Vsp+Vhe+Vlu+Vbo+Vbr+Vmu+Vad+VguWall+Vbl)

        CO = C_OUTPUT*60
        Qad = FQad*CO
        Qbo = FQbo*CO
        Qbr = FQbr*CO
        Qgu = FQgu*CO
        Qhe = FQhe*CO
        Qki = FQki*CO
        Qli = FQli*CO
        Qmu = FQmu*CO
        Qsp = FQsp*CO
        Qha = Qli - (Qgu+Qsp)
        Qtot = Qli+Qki+Qbo+Qhe+Qmu+Qad+Qbr
        Qre = CO - Qtot
        Qlu = CO

        SA_abs = pi*L*d*PF*VF*MF*1e-4
        SA_basal = pi*L*d*PF*VF*1e-4
        MA = 10^logP
        MW_eff = MW - (3*17)
        Peff = fperm*A*(((MW_eff^(-alpha-beta))*MA)/((MW_eff^(-alpha)) + B*(MW_eff^(-beta))*MA) * 1e-2 * 3600)
        kd = fdis*Peff*SA_abs*1000/VguLumen
        ka = fabs*Peff*SA_basal*1000/VguWall
        kt = 1/ITT

        scale_factor_H = MPPGL*Vli*1000
        scale_factor_G = MPPGI*VguWall*1000

        CLintHep = ((VmaxH/KmH)*scale_factor_H*60*1e-6)/fumic
        CLintGut = ((VmaxG/KmG)*scale_factor_G*60*1e-6)/fumic
        #CLintHep = CLintHep/fumic
        #CLintGut = CLintGut/fumic

        CLrenal = 0.096

        f = 1
    end

    @dynamics begin
        GUTLUMEN' = -kd*VguLumen*(f*(GUTLUMEN/VguLumen) + (1-f)*S_lumen) -
            kt*GUTLUMEN
        GUTWALL' = kd*VguLumen*(f*(GUTLUMEN/VguLumen) + (1-f)*S_lumen) -
            ka*GUTWALL - CLintGut*(GUTWALL/VguWall)
        GUT' = ka*GUTWALL + Qgu*((ART/Var) - (GUT/VguWall)/(Kpgu/BP))
        ADIPOSE' = Qad*((ART/Var) - (ADIPOSE/Vad)/(Kpad/BP))
        BRAIN' = Qbr*((ART/Var) - (BRAIN/Vbr)/(Kpbr/BP))
        HEART' = Qhe*((ART/Var) - (HEART/Vhe)/(Kphe/BP))
        KIDNEY' = Qki*((ART/Var) - (KIDNEY/Vki)/(Kpki/BP)) -
            CLrenal*(((KIDNEY/Vki)*fup)/(Kpki/BP))
        LIVER' = Qgu*((GUT/VguWall)/(Kpgu/BP)) + Qsp*((SPLEEN/Vsp)/(Kpsp/BP)) +
            Qha*(ART/Var) - Qli*((LIVER/Vli)/(Kpli/BP)) -
            CLintHep*(((LIVER/Vli)*fup)/(Kpli/BP))
        LUNG' = Qlu*((VEN/Vve) - (LUNG/Vlu)/(Kplu/BP))
        MUSCLE' = Qmu*((ART/Var) - (MUSCLE/Vmu)/(Kpmu/BP))
        SPLEEN' = Qsp*((ART/Var) - (SPLEEN/Vsp)/(Kpsp/BP))
        BONE' = Qbo*((ART/Var) - (BONE/Vbo)/(Kpbo/BP))
        REST' = Qre*((ART/Var) - (REST/Vre)/(Kpre/BP))
        VEN' = Qad*((ADIPOSE/Vad)/(Kpad/BP)) + Qbr*((BRAIN/Vbr)/(Kpbr/BP)) +
            Qhe*((HEART/Vhe)/(Kphe/BP)) + Qki*((KIDNEY/Vki)/(Kpki/BP)) +
            Qli*((LIVER/Vli)/(Kpli/BP)) + Qmu*((MUSCLE/Vmu)/(Kpmu/BP)) +
            Qbo*((BONE/Vbo)/(Kpbo/BP)) + Qre*((REST/Vre)/(Kpre/BP)) -
            Qlu*(VEN/Vve)
        ART' = Qlu*((LUNG/Vlu)/(Kplu/BP) - (ART/Var))
    end

    @derived begin
        gutlumen = GUTLUMEN
        gutwall = GUTWALL
        gut = GUT
        adipose = ADIPOSE
        brain = BRAIN
        heart = HEART
        bone = BONE
        kidney = KIDNEY
        liver = LIVER
        lung = LUNG
        muscle = MUSCLE
        spleen = SPLEEN
        rest = REST
        art = ART
        ven = VEN
        Cvenn = VEN/Vve
        nca := @nca Cvenn
        auc =  NCA.auc(nca, auctype = :last)
        cmax = NCA.cmax(nca)
    end
end

#parameters - tissue:plasma partition coefficients based on: Poulin and Theil
p = (fup = 0.42, fumic = 0.711, WEIGHT = 73, MPPGL = 30.3, MPPGI = 0,
    C_OUTPUT = 6.5, VmaxH = 40, VmaxG = 40, KmH = 9.3, KmG = 9.3, BP = 1,
    Kpad = 9.89, Kpbo = 7.91, Kpbr = 7.35, Kpgu = 5.82, Kphe = 1.95, Kpki = 2.9,
    Kpli = 4.66, Kplu = 0.83, Kpmu = 2.94, Kpsp = 2.96, Kpre = 4, MW = 349.317,
    logP = 2.56, S_lumen = 0.39*1000, L = 280, d = 2.5, PF = 1.57, VF = 6.5,
    MF = 13, ITT = 3.32, A = 7440, B = 1e7, alpha = 0.6, beta = 4.395, fabs = 1,
    fdis = 1, fperm = 1, Vad = 18.2, Vbo = 10.5, Vbr = 1.45, VguWall = 0.65,
    VguLumen = 0.35, Vhe = 0.33, Vki = 0.31, Vli = 1.8, Vlu = 0.5, Vmu = 29,
    Vsp = 0.15, Vbl = 5.6, FQad = 0.05, FQbo = 0.05, FQbr = 0.12, FQgu = 0.16,
    FQhe = 0.04, FQki = 0.19, FQli = 0.255, FQmu = 0.17, FQsp = 0.03)

#Weight = 73 kg, regimen = 4mg/kg q12h IV
regimen = DosageRegimen(292, time = 0, addl=13, ii=12, cmt=14, rate = 292, ss = 1)

sub1 = Subject(id=1,evs=regimen)

#predictions based on Poulin and Theil http://jpharmsci.org/article/S0022-3549(16)30889-9/fulltext
pred_PT = simobs(model, sub1, p, obstimes=0:0.1:12)

#parameters - tissue:plasma partition coefficients based on: Berezhkowskiy, Leonid M (2004)
p_Berez = (p..., Kpad = 10.7, Kpbo = 7.33, Kpbr = 6.87, Kpgu = 5.82, Kphe = 1.92,
        Kpki = 2.8, Kpli = 4.4, Kplu = 0.91, Kpmu = 2.83, Kpsp = 2.85, Kpre = 4) #Kpgu ??

#predictions based on Berez
pred_Berez = simobs(model, sub1, p_Berez, obstimes=0:0.1:12) #slightly different numbers - could be due to Kpgu?

#parameters - tissue:plasma partition coefficients based on: Rodgers and Rowland http://jpharmsci.org/article/S0022-3549(16)31789-0/fulltext and http://jpharmsci.org/article/S0022-3549(16)32034-2/fulltext
p_RR = (p..., Kpad = 11.5, Kpbo = 2.86, Kpbr = 6.34, Kpgu = 6.69, Kphe = 2.97,
        Kpki = 3.26, Kpli = 3.55, Kplu = 4.27, Kpmu = 2.17, Kpsp = 2.02, Kpre = 4)

#predictions based on Rodgers and Rowland
pred_RR = simobs(model, sub1, p_RR, obstimes = 0:0.1:12) #slightly different numbers

#parameters - tissue:plasma partition coefficients based on: Schmitt https://www.sciencedirect.com/science/article/pii/S0887233307002573?via%3Dihub
p_Schmitt = (p..., Kpad = 141, Kpbo = 3.88, Kpbr = 17.4, Kpgu = 5.82, Kphe = 17.7,
        Kpki = 10.2, Kpli = 13.1, Kplu = 6.79, Kpmu = 2.5, Kpsp = 4.12, Kpre = 4)

#predictions based on Schmitt
pred_Schmitt = simobs(model, sub1, p_Schmitt, obstimes = 0:0.1:12) #slightly different numbers

#parameters - tissue:plasma partition coefficients based on: PK-Sim https://www.tandfonline.com/doi/abs/10.1517/17425255.1.1.159
p_pksim = (p..., Kpad = 122, Kpbo = 42, Kpbr = 17.4, Kpgu = 5.82, Kphe = 16.1,
        Kpki = 8.84, Kpli = 11.5, Kplu = 2.49, Kpmu = 2.93, Kpsp = 3.43, Kpre = 4)

#predictions based on PK-Sim
pred_pksim = simobs(model, sub1, p_pksim, obstimes = 0:0.1:12) #slightly different numbers


#csv import
#observed (digitized) data from fig3a in the ZT paper
obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/obs.csv"))
#ZT predicted data (digitized) from fig3a in the ZT paper
ZT = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/ZT.csv"))


#plot:
#ZT predicted vs observed vs PT/Berez/RR/Schmitt/PKSim predicted for Adult 4mg/kg IV
plot(pred_PT.times, pred_PT[:Cvenn], linecolor = :black, label = "PT", title = "Adult 4 mg/kg IV")
plot!(pred_Berez.times, pred_Berez[:Cvenn], line = [:dash], linecolor = :black, label = "Berez")
plot!(pred_RR.times, pred_RR[:Cvenn], line = [:dot], linecolor = :black, label = "RR")
plot!(pred_Schmitt.times, pred_Schmitt[:Cvenn], line = :dashdot, linecolor = :black, label = "Schmitt")
plot!(pred_pksim.times, pred_pksim[:Cvenn], line = :dashdotdot, linecolor = :black, label = "PK Sim")
plot!(obs.time, obs[:obs], seriestype=:scatter, yerror = obs[:sd], markercolor = :black, seriescolor = :black, label = "observed")
plot!(ZT.time, ZT[:ZT], linecolor = :lightgray, label = "ZT")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")

png("adult_1.png")
#Table 1: This chunk tables out the prediction errors of different calculation methods
#pad observed data till 12 hours

#filter observed data for time > 7
df_temp = @from i in obs begin
    @where i.time > 7
    @select {i.time, i.obs, i.sd}
    @collect DataFrame
end

mod = coef(lm(@formula(obs ~ time), df_temp))

#adding time 12 to data frame and predicting using lm
st = DataFrame(time = [12])
st2 = GLM.predict((lm(@formula(obs ~ time), df_temp)), st)

df_temp = DataFrame(Column1 = 14, time = 12, obs = st2, sd = 0)
append!(obs, df_temp)

## calculate AUCs for all methods
auc_obs = NCA.auc(obs[:obs], obs[:time], interval = (0, 12))
auc_PT = NCA.auc(pred_PT[:Cvenn], pred_PT.times, interval = (0,12))
auc_Berez = NCA.auc(pred_Berez[:Cvenn], pred_Berez.times, interval = (0,12))
auc_RR = NCA.auc(pred_RR[:Cvenn], pred_RR.times, interval = (0,12))
auc_Schmitt = NCA.auc(pred_Schmitt[:Cvenn], pred_Schmitt.times, interval = (0,12))
auc_pksim = NCA.auc(pred_pksim[:Cvenn], pred_pksim.times, interval = (0,12))

## calculate cmax for all methods
cmax_obs = NCA.cmax(obs[:obs], obs[:time])
cmax_PT = NCA.cmax(pred_PT[:Cvenn], pred_PT.times)
cmax_Berez = NCA.cmax(pred_Berez[:Cvenn], pred_Berez.times)
cmax_RR = NCA.cmax(pred_RR[:Cvenn], pred_RR.times)
cmax_Schmitt = NCA.cmax(pred_Schmitt[:Cvenn], pred_Schmitt.times)
cmax_pksim = NCA.cmax(pred_pksim[:Cvenn], pred_pksim.times)

#need to find out how to mutate


## Figure 4a; Model 1 with 4 mg/kg IV infusion - same as before
regimen2 = DosageRegimen(292, time = 0, addl=7, ii=12, cmt=14, rate = 292, ss = 1)

sub2 = Subject(id=1,evs=regimen2)
sim = simobs(model, sub2, p_RR, obstimes = 0:0.1:12) #slightly different numbers

#csv import attempt
obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/obs.csv"))
ZT = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/ZT.csv"))

plot(sim.times, sim[:Cvenn], linecolor = :black, label = "model", title = "a Adult 4 mg/kg IV")
plot!(obs.time, obs[:obs], seriestype=:scatter, yerror = obs[:sd], markercolor = :black, seriescolor = :black, label = "observed")
plot!(ZT.time, ZT[:ZT], linecolor = :blue, label = "ZT")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")

png("adult_1.png")

## Figure 3b; Model 2 with 4 mg/kg IV infusion
obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4a_obs.csv")) #obs data from fig 3a in the ZT paper
ZT = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4a_ZT.csv")) #ZT predictions from fig 3a in the ZT paper

#Model 2 - pediatrics
fixeffs_ped = (fup = 0.42, fumic = 0.711, WEIGHT = 19, MPPGL = 26, MPPGI = 0,
    C_OUTPUT = 3.4, VmaxH = 120.5, VmaxG = 120.5, KmH = 11, KmG = 11, BP = 1,
    Kpad = 9.89, Kpbo = 7.91, Kpbr = 7.35, Kpgu = 5.82, Kphe = 1.95, Kpki = 2.9,
    Kpli = 4.66, Kplu = 0.83, Kpmu = 2.94, Kpsp = 2.96, Kpre = 4, MW = 349.317,
    logP = 2.56, S_lumen = 0.39*1000, L = 170, d = 2.5, PF = 1.57, VF = 6.5,
    MF = 13, ITT = 3.32, A = 7440, B = 1e7, alpha = 0.6, beta = 4.395, fabs = 1,
    fdis = 1, fperm = 1, Vad = 5.5, Vbo = 2.43, Vbr = 1.31, VguWall = 0.22,
    VguLumen = 0.117, Vhe = 0.085, Vki = 0.11, Vli = 0.467, Vlu = 0.125, Vmu = 5.6,
    Vsp = 0.05, Vbl = 1.5, FQad = 0.05, FQbo = 0.05, FQbr = 0.12, FQgu = 0.16,
    FQhe = 0.04, FQki = 0.19, FQli = 0.255, FQmu = 0.17, FQsp = 0.03)

#updating pediatric model with RR Kps
p_RR_ped = (fixeffs_ped..., Kpad = 11.5, Kpbo = 2.86, Kpbr = 6.34, Kpgu = 6.69, Kphe = 2.97,
        Kpki = 3.26, Kpli = 3.55, Kplu = 4.27, Kpmu = 2.17, Kpsp = 2.02, Kpre = 4)

#Pediatric BW = 19 kg, dose = 4 mg/kg, rate = 3 mg/kg, route = IV
regimen3 = DosageRegimen(76, time = 0, addl=7, ii=12, cmt=14, rate = 57, ss = 1)
sub3 = Subject(id=1,evs=regimen3)
sim = simobs(model, sub3, p_RR_ped, obstimes = 0:0.1:12)

#Pediatric 4 mg/kg IV plot
plot(sim.times, sim[:Cvenn], linecolor = :black, label = "model", title = "b Pediatric 4 mg/kg IV")
plot!(obs.time, obs[:obs], seriestype=:scatter, yerror = obs[:sd], markercolor = :black, seriescolor = :black, label = "observed")
plot!(ZT.time, ZT[:ZT], linecolor = :blue, label = "ZT")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")

png("ped_1.png")
#Figure 3c; Model 1 with 200 mg PO
obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig3b_obs.csv")) #observed data from fig 3b in the ZT paper
ZT = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig3b_ZT.csv")) #predicted data from fig 3b in the ZT paper

#200 mg PO dose
regimen4 = DosageRegimen(200, time = 0, addl=7, ii=12, cmt=1, ss = 1)
sub4 = Subject(id=1, evs=regimen4)
sim = simobs(model, sub4, p_RR, obstimes = 0:0.1:12)

plot(sim.times, sim[:Cvenn], linecolor = :black, label = "model 1", title = "c Adult 200 mg PO")
plot!(obs.time, obs[:obs], seriestype=:scatter, yerror = obs[:sd], markercolor = :black, seriescolor = :black, label = "observed")
plot!(ZT.time, ZT[:ZT], linecolor = :blue, label = "ZT")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")

png("adult_1_po.png")

#Figure 3d; Model 2 with 4 mg/kg PO
obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4d_obs.csv")) #observed data from fig 4d in ZT paper
ZT = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4d_ZT.csv")) #predicted data from fig 4d in ZT paper

#Weight = 19 kg, Dose = 4 mg/kg PO
regimen5 = DosageRegimen(76, time = 0, addl=7, ii=12, cmt=1, ss = 1)
sub5 = Subject(id=1, evs=regimen5)
sim = simobs(model, sub5, p_RR_ped, obstimes = 0:0.1:12)

plot(sim.times, sim[:Cvenn], linecolor = :black, label = "model 1", title = "d Pediatric 4 mg/kg PO")
plot!(obs.time, obs[:obs], seriestype=:scatter, yerror = obs[:sd], markercolor = :black, seriescolor = :black, label = "observed")
plot!(ZT.time, ZT[:ZT], linecolor = :blue, label = "ZT")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")

png("ped_1_po.png")

### sensitivity analyses for absorption parameters: permeability, intestinal transit time and solubility ###
#####SENSITIVITY ANALYSIS STARTS HERE
fixeffs = p

regimen_s = DosageRegimen(200, time = 0, addl=13, ii=12, cmt=1, ss = 1)

sub_s = Subject(id=1,evs=regimen_s)

function generate_population(events,nsubs=4)
  pop = Population(map(i -> Subject(id=i,evs=events),1:nsubs))
  return pop
end

sub_p = generate_population(regimen_s)

#Kpmu sensitivity analysis for abstract
_Kpmu1 = (fixeffs..., Kpmu = 1.5)
_Kpmu2 = (fixeffs..., Kpmu = 3)
_Kpmu3 = (fixeffs..., Kpmu = 6)

pred_kpmu1 = simobs(model, sub_s, _Kpmu1)
pred_kpmu2 = simobs(model, sub_s, _Kpmu2)
pred_kpmu3 = simobs(model, sub_s, _Kpmu3)

p1 = plot(pred_kpmu1.times, pred_kpmu1[:Cvenn], linecolor = :black, label = "1.5", legendtitle = "Kpmu", legendfontsize = 12)
plot!(pred_kpmu2.times, pred_kpmu2[:Cvenn], linecolor = :black, line = :dash, label = "3")
plot!(pred_kpmu3.times, pred_kpmu3[:Cvenn], linecolor = :black, line = :dot, label = "6")
xlabel!("Time (h)")
ylabel!("Plasma concentration (mg/L)")
p1

#png("kpmu.png")


#BP sensitivity analysis for abstract
_BP1 = (fixeffs..., BP = 0.5)
_BP2 = (fixeffs..., BP = 1)
_BP3 = (fixeffs..., BP = 2)

pred_BP1 = simobs(model, sub_s, _BP1)
pred_BP2 = simobs(model, sub_s, _BP2, obstimes=0:0.1:12)
pred_BP3 = simobs(model, sub_s, _BP3, obstimes=0:0.1:12)

p2 = plot(pred_BP1.times, pred_BP1[:Cvenn], linecolor = :black, label = "0.5", legendtitle = "BP", legendfontsize = 12)
plot!(pred_BP2.times, pred_BP2[:Cvenn], linecolor = :black, line = :dash, label = "1")
plot!(pred_BP3.times, pred_BP3[:Cvenn], linecolor = :black, line = :dot, label = "2")
xlabel!("Time (h)")
ylabel!("Plasma concentration (mg/L)")

#png("bp.png")

#parameters - tissue:plasma partition coefficients based on: Poulin and Theil
p = (fup = 0.42, fumic = 0.711, WEIGHT = 73, MPPGL = 30.3, MPPGI = 0,
    C_OUTPUT = 6.5, VmaxH = 40, VmaxG = 40, KmH = 9.3, KmG = 9.3, BP = 1,
    Kpad = 9.89, Kpbo = 7.91, Kpbr = 7.35, Kpgu = 5.82, Kphe = 1.95, Kpki = 2.9,
    Kpli = 4.66, Kplu = 0.83, Kpmu = 2.94, Kpsp = 2.96, Kpre = 4, MW = 349.317,
    logP = 2.56, S_lumen = 0.39*1000, L = 280, d = 2.5, PF = 1.57, VF = 6.5,
    MF = 13, ITT = 3.32, A = 7440, B = 1e7, alpha = 0.6, beta = 4.395, fabs = 1,
    fdis = 1, fperm = 1, Vad = 18.2, Vbo = 10.5, Vbr = 1.45, VguWall = 0.65,
    VguLumen = 0.35, Vhe = 0.33, Vki = 0.31, Vli = 1.8, Vlu = 0.5, Vmu = 29,
    Vsp = 0.15, Vbl = 5.6, FQad = 0.05, FQbo = 0.05, FQbr = 0.12, FQgu = 0.16,
    FQhe = 0.04, FQki = 0.19, FQli = 0.255, FQmu = 0.17, FQsp = 0.03)

sobol = gsa(model,
            sub_p,
            fixeffs,
            DiffEqSensitivity.Sobol(N=2000,order=[0,1], nboot = 50),
            [:auc],
           (Kpmu = 1.5, BP=0.5, S_lumen=0.39*500, fperm=0.5, MPPGI=1.212/2, ITT=3.32/2),
           (Kpmu = 6, BP=2, S_lumen=0.39*2000, fperm=2, MPPGI=1.212*2, ITT=3.32*2))

serialize("sobolresult.jl", sobol)

#changing S_lumen
_fixeffs = (fixeffs...,S_lumen = 0.39*500)
_fixeffs2 = (fixeffs...,S_lumen = 0.39*1000)
_fixeffs3 = (fixeffs...,S_lumen = 0.39*2000)

pred_1 = simobs(model, sub_p, _fixeffs)
pred_2 = simobs(model, sub_s, _fixeffs2, obstimes=0:0.1:12)
pred_3 = simobs(model, sub_s, _fixeffs3, obstimes=0:0.1:12)

plot(pred_1.times, pred_1[:Cvenn], linecolor = :black, label = "0.195", title = "Sint (mg/mL)")
plot!(pred_2.times, pred_2[:Cvenn], linecolor = :black, line = :dash, label = "0.39")
plot!(pred_3.times, pred_3[:Cvenn], linecolor = :black, line = :dot, label = "0.78")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")

#changing fperm
_fixeffs4 = (fixeffs...,fperm = 0.5)
_fixeffs5 = (fixeffs...,fperm = 1)
_fixeffs6 = (fixeffs...,fperm = 2)


pred_4 = simobs(model, sub_s, _fixeffs4)
pred_5 = simobs(model, sub_s, _fixeffs5, obstimes=0:0.1:12)
pred_6 = simobs(model, sub_s, _fixeffs6, obstimes=0:0.1:12)

#my_nca = read_nca(DataFrame(pred_4), id = :id, conc = :Cvenn, time = :time, amt = :amt)


plot(pred_4.times, pred_4[:Cvenn], linecolor = :black, label = "0.073", title = "Peff x 10-4 (cm/s)")
plot!(pred_5.times, pred_5[:Cvenn], linecolor = :black, line = :dash, label = "0.145")
plot!(pred_6.times, pred_6[:Cvenn], linecolor = :black, line = :dot, label = "0.29")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")

png("peff.png")


#changing ITT
_fixeffs7 = (fixeffs...,ITT = 3.32/2)
_fixeffs8 = (fixeffs...,ITT = 3.32)
_fixeffs9 = (fixeffs...,ITT = 3.32*2)

pred_7 = simobs(model, sub_s, _fixeffs7, obstimes=0:0.1:12)
pred_8 = simobs(model, sub_s, _fixeffs8, obstimes=0:0.1:12)
pred_9 = simobs(model, sub_s, _fixeffs9, obstimes=0:0.1:12)

plot(pred_7.times, pred_7[:Cvenn], linecolor = :black, label = "1.66", title = "ITT (h)")
plot!(pred_8.times, pred_8[:Cvenn], linecolor = :black, line = :dash, label = "3.32")
plot!(pred_9.times, pred_9[:Cvenn], linecolor = :black, line = :dot, label = "6.64")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")

#changing MPPGI
_fixeffs10 = (fixeffs...,MPPGI = 1.212/2)
_fixeffs11 = (fixeffs...,MPPGI = 1.212)
_fixeffs12 = (fixeffs...,MPPGI = 1.212*2)

pred_10 = simobs(model, sub_s, _fixeffs10, obstimes=0:0.1:12)
pred_11 = simobs(model, sub_s, _fixeffs11, obstimes=0:0.1:12)
pred_12 = simobs(model, sub_s, _fixeffs12, obstimes=0:0.1:12)

plot(pred_10.times, pred_10[:Cvenn], linecolor = :black, label = "0.035", title = "CLgu (mL/min/kg)")
plot!(pred_11.times, pred_11[:Cvenn], linecolor = :black, line = :dash, label = "0.07")
plot!(pred_12.times, pred_12[:Cvenn], linecolor = :black, line = :dot, label = "0.14")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")

png("clgu.png")


###Figure 7 - reproduces figure 5 plots
## calculate fperm (permeability factor)
VguWall = 0.65
VguLumen = 0.35
MW = 349.317  #(g/mol)
logP = 2.56  #log10 octanol oil:water partition coefficient; will be used as proxy for membrane affinity; preferably we will have phospholipid bilayer:water partition instead
S_lumen = 0.39*1000  #(mg/L) voriconazole intestinal lumen solubility https://www.ncbi.nlm.nih.gov/pubmed/24557773
L = 280  #(cm) small intestine length; from ICRP Publication 89
d = 2.5  #(cm) diameter of small intestine lumen
PF = 1.57  #3  //2.29  //average of 1.57 and 3; 1.57  //plicae circulare factor https://www.ncbi.nlm.nih.gov/pubmed/24694282
VF = 6.5  #villi factor
MF = 13  #microvilli factor
A = 7440  #this and the rest of parameters are constants in the permeability calculation equation https://www.ncbi.nlm.nih.gov/pubmed/15267240
B = 1e7
alpha = 0.6
beta = 4.395
fabs = 1
fperm = 1
SA_basal = pi*L*d*PF*VF*1e-4
MA = 10^logP
MW_eff = MW - (3*17)
Peff = fperm*A*((MW_eff^(-alpha-beta))*MA/(MW_eff^(-alpha) + B*(MW_eff^(-beta))*MA) * 1e-2 * 3600)
ka = fabs*Peff*SA_basal*1000/VguWall  #(h-1)
fperm = 0.849/ka  #0.849 is the reported absorption rate constant

p_RR = (fixeffs..., Kpad = 11.5, Kpbo = 2.86, Kpbr = 6.34, Kpgu = 6.69, Kphe = 2.97,
        Kpki = 3.26, Kpli = 3.55, Kplu = 4.27, Kpmu = 2.17, Kpsp = 2.02, Kpre = 4)

model3 = (p_RR..., MPPGI=30.3/25)
model4 = (p_RR_ped..., MPPGI = 26/25)
model5 = (p_RR..., fperm = fperm)
model6 = (p_RR_ped..., fperm = fperm)
model7 = (p_RR..., MPPGI = 30.3/25, fperm = fperm)
model8 = (p_RR_ped..., MPPGI = 26/25, fperm = fperm)


#Figure 5a; Model 3 with 200 mg PO
obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig3b_obs.csv")) #observed data from fig 3b in ZT paper
ZT = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig3b_ZT.csv")) #predicted data from fig 3b in ZT paper

#200 mg PO
regimen = DosageRegimen(200, time = 0, addl=7, ii=12, cmt=1, ss = 1)

sub = Subject(id=1,evs=regimen)

sim = simobs(model, sub, p_RR, obstimes=0:0.1:12)
sim_mppgi = simobs(model, sub, model3, obstimes=0:0.1:12) #MPPGI = 30.3/25
sim_fperm = simobs(model, sub, model5, obstimes=0:0.1:12) #calculated fperm
sim_mppgi_fperm = simobs(model, sub, model7, obstimes=0:0.1:12) #MPPGI = 30.3/25 & calculated fperm

plot(sim.times, sim[:Cvenn], line = [:dash],linecolor = :black, label = "model 1", title = "Adult 200mg PO")
plot!(sim_mppgi.times, sim_mppgi[:Cvenn], line = [:dashdot], linecolor = :black, label = "model 3")
plot!(sim_fperm.times, sim_fperm[:Cvenn], line = [:dot], linecolor = :black, label = "model 5")
plot!(sim_mppgi_fperm.times, sim_mppgi_fperm[:Cvenn], linecolor = :black, label = "model 7")
plot!(obs.time, obs[:obs], seriestype=:scatter, yerror = obs[:sd], markercolor = :black, seriescolor = :black, label = "observed")
plot!(ZT.time, ZT[:ZT], linecolor = :lightgray, label = "ZT")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")

plot(sim.times, sim[:Cvenn], line = [:dash],linecolor = :black, label = "initial model", title = "Adult 200mg PO")
plot!(sim_mppgi_fperm.times, sim_mppgi_fperm[:Cvenn], linecolor = :black, label = "model w/Clgu")
plot!(obs.time, obs[:obs], seriestype=:scatter, yerror = obs[:sd], markercolor = :black, seriescolor = :black, label = "observed")
plot!(ZT.time, ZT[:ZT], linecolor = :blue, label = "ZT")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")

png("adult_po_clgu.png")

#Figure 5b; Model 4 with 4 mg/kg IV infusion
obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4d_obs.csv")) #observed data from fig 4d in ZT paper
ZT_Gu = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4e_ZT.csv")) #predicted data from fig 4d in ZT paper
ZT = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4d_ZT.csv")) #predicted data from fig 4d in the ZT paper

#Weight = 19 kg, DOse = 4 mg/kg PO - pediatric
regimen = DosageRegimen(76, time = 0, addl=7, ii=12, cmt=1, ss = 1)
sub = Subject(id=1,evs=regimen)

sim = simobs(model, sub, p_RR_ped, obstimes=0:0.1:12)
sim_mppgi = simobs(model, sub, model4, obstimes=0:0.1:12)
sim_fperm = simobs(model, sub, model6, obstimes=0:0.1:12)
sim_mppgi_fperm = simobs(model, sub, model8, obstimes=0:0.1:12)


plot(sim.times, sim[:Cvenn], line = [:dash],linecolor = :black, label = "model 1", title = "Pediatric 4mg/kg PO")
plot!(sim_mppgi.times, sim_mppgi[:Cvenn], line = [:dashdot], linecolor = :black, label = "model 4")
plot!(sim_fperm.times, sim_fperm[:Cvenn], line = [:dot], linecolor = :black, label = "model 6")
plot!(sim_mppgi_fperm.times, sim_mppgi_fperm[:Cvenn], linecolor = :black, label = "model 8")
plot!(obs.time, obs[:obs], seriestype=:scatter, yerror = obs[:sd], markercolor = :black, seriescolor = :black, label = "observed")
plot!(ZT.time, ZT[:ZT], linecolor = :lightgray, label = "ZT")
plot!(ZT_Gu.time, ZT_Gu[:ZT], linecolor = :lightgray, line = :dash, label ="ZT_Gu")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")

plot(sim.times, sim[:Cvenn], line = [:dash],linecolor = :black, label = "initial model", title = "Pediatric 4mg/kg PO")
plot!(sim_mppgi_fperm.times, sim_mppgi_fperm[:Cvenn], linecolor = :black, label = "model w/Clgu")
plot!(obs.time, obs[:obs], seriestype=:scatter, yerror = obs[:sd], markercolor = :black, seriescolor = :black, label = "observed")
plot!(ZT.time, ZT[:ZT], linecolor = :blue, label = "ZT")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")

png("ped_po_clgu.png")

###Table 2
obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/obs.csv"))
ZT = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/ZT.csv"))

df_temp = @from i in obs begin
    @where i.time > 7
    @select {i.time, i.obs, i.sd}
    @collect DataFrame
end

mod = coef(lm(@formula(obs ~ time), df_temp))

st = DataFrame(time = [12])
st2 = predict((lm(@formula(obs ~ time), df_temp)), st)

df_temp = DataFrame(Column1 = 14, time = 12, obs = st2, sd = 0)
append!(obs, df_temp)

#Weight = 73 kg, Dose = 4 mg/kg IV
regimen = DosageRegimen(292, time = 0, addl=13, ii=12, rate = 292, cmt=14, ss = 1)
sub = Subject(id=1,evs=regimen)

pred = simobs(model, sub, p_RR, obstimes=0:0.1:12)

auc_obs = NCA.auc(obs[:obs], obs[:time], interval = (0, 12))
auc_ZT = NCA.auc(ZT[:ZT], ZT[:time], interval = (0, 12))
auc_pred = NCA.auc(pred[:Cvenn], pred.times, interval = (0,12))

cmax_obs = NCA.cmax(obs[:obs], obs[:time])
cmax_ZT = NCA.cmax(ZT[:ZT], ZT[:time])
cmax_pred = NCA.cmax(pred[:Cvenn], pred.times)

#need to make into dataframe and then use mutate

#pediatric load data
obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4a_obs.csv"))
ZT = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4a_ZT.csv"))

#Weight = 19 kg, Dose = 4 mg/kg IV
regimen = DosageRegimen(76, time = 0, addl=13, ii=12, cmt = 14, rate = 76, ss = 1)
sub = Subject(id=1, evs=regimen)

pred = simobs(model, sub, p_RR_ped, obstimes=0:0.1:12)

auc_obs = NCA.auc(obs[:obs], obs[:time], interval = (0, 12))
auc_ZT = NCA.auc(ZT[:ZT], ZT[:time], interval = (0, 12))
auc_pred = NCA.auc(pred[:Cvenn], pred.times, interval = (0,12))

cmax_obs = NCA.cmax(obs[:obs], obs[:time])
cmax_ZT = NCA.cmax(ZT[:ZT], ZT[:time])
cmax_pred = NCA.cmax(pred[:Cvenn], pred.times)

#need to make into dataframe and then use mutate

#Table 3
#adult load data
obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig3b_obs.csv"))
ZT = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig3b_ZT.csv"))

df_temp = @from i in obs begin
    @where i.time > 7
    @select {i.time, i.obs, i.sd}
    @collect DataFrame
end

mod = coef(lm(@formula(obs ~ time), df_temp))

st = DataFrame(time = [12])
st2 = predict((lm(@formula(obs ~ time), df_temp)), st)

df_temp = DataFrame(Column1 = 14, time = 12, obs = st2, sd = 0)
append!(obs, df_temp)

#Dose = 200 mg PO
regimen = DosageRegimen(200, time = 0, addl=13, ii=12, cmt = 1, ss = 1)
sub = Subject(id=1, evs=regimen)

pred_init = simobs(model, sub, p_RR, obstimes=0:0.1:12)
pred_mppgi = simobs(model, sub, model3, obstimes=0:0.1:12)
pred_fperm = simobs(model, sub, model5, obstimes=0:0.1:12)
pred_mppgi_fperm = simobs(model, sub, model7, obstimes=0:0.1:12)

auc_obs = NCA.auc(obs[:obs], obs[:time], interval = (0, 12))
auc_ZT = NCA.auc(ZT[:ZT], ZT[:time], interval = (0, 12))
auc_pred_init = NCA.auc(pred_init[:Cvenn], pred_init.times, interval = (0,12))
auc_pred_mppgi = NCA.auc(pred_mppgi[:Cvenn], pred_mppgi.times, interval = (0,12))
auc_pred_fperm = NCA.auc(pred_fperm[:Cvenn], pred_fperm.times, interval = (0,12))
auc_pred_mppgi_fperm = NCA.auc(pred_mppgi_fperm[:Cvenn], pred_mppgi_fperm.times, interval = (0,12))


cmax_obs = NCA.cmax(obs[:obs], obs[:time])
cmax_ZT = NCA.cmax(ZT[:ZT], ZT[:time])
cmax_pred_init = NCA.cmax(pred_init[:Cvenn], pred_init.times)
cmax_pred_mppgi = NCA.cmax(pred_mppgi[:Cvenn], pred_mppgi.times)
cmax_pred_fperm = NCA.cmax(pred_fperm[:Cvenn], pred_fperm.times)
cmax_pred_mppgi_fperm = NCA.cmax(pred_mppgi_fperm[:Cvenn], pred_mppgi_fperm.times)

#need to make into dataframe and then use mutate

#load pediatric data
obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4d_obs.csv"))
ZT = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4e_ZT.csv"))
df_temp = DataFrame(Column1 = 0, time = 0, obs = 0, sd = 0)
append!(obs, df_temp)
sort!(obs)


#Weight = 19 kg, Dose = 4mg/kg PO
regimen = DosageRegimen(76, time = 0, addl=13, ii=12, cmt = 1, ss = 1)
sub = Subject(id=1, evs=regimen)

pred_init = simobs(model, sub, p_RR_ped, obstimes=0:0.1:12)
pred_mppgi = simobs(model, sub, model4, obstimes=0:0.1:12)
pred_fperm = simobs(model, sub, model6, obstimes=0:0.1:12)
pred_mppgi_fperm = simobs(model, sub, model8, obstimes=0:0.1:12)

auc_obs = NCA.auc(obs[:obs], obs[:time], interval = (0, 12)) #auc value different by 0.7
auc_ZT = NCA.auc(ZT[:ZT], ZT[:time], interval = (0, 12))
auc_pred_init = NCA.auc(pred_init[:Cvenn], pred_init.times, interval = (0,12))
auc_pred_mppgi = NCA.auc(pred_mppgi[:Cvenn], pred_mppgi.times, interval = (0,12))
auc_pred_fperm = NCA.auc(pred_fperm[:Cvenn], pred_fperm.times, interval = (0,12))
auc_pred_mppgi_fperm = NCA.auc(pred_mppgi_fperm[:Cvenn], pred_mppgi_fperm.times, interval = (0,12))

obs

cmax_obs = NCA.cmax(obs[:obs], obs[:time])
cmax_ZT = NCA.cmax(ZT[:ZT], ZT[:time])
cmax_pred_init = NCA.cmax(pred_init[:Cvenn], pred_init.times)
cmax_pred_mppgi = NCA.cmax(pred_mppgi[:Cvenn], pred_mppgi.times)
cmax_pred_fperm = NCA.cmax(pred_fperm[:Cvenn], pred_fperm.times)
cmax_pred_mppgi_fperm = NCA.cmax(pred_mppgi_fperm[:Cvenn], pred_mppgi_fperm.times)

#need to make into dataframe and then use mutate


###IGNORE code below
#Figure 8

#adults
obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig3b_obs.csv"))

Random.seed!(23142)
nSampl = 1000
x = rand(Normal(log(30.3/25),0.33),nSampl)
mppgiSampl = exp.(x)
fpermSampl = exp.(rand(Normal(log(fperm), 0.4), nSampl))  #lognormal; CV = 40%
lpars = DataFrame(MPPGI = mppgiSampl, fperm = fpermSampl)
lpars = hcat(mppgiSampl, fpermSampl)

function Foo(x,y)
    mod = (p_RR...,MPPGI = x, fperm = y)
    regimen = DosageRegimen(200, time = 0, addl=13, ii=12, cmt = 1, ss = 1)
    sub = Subject(id=1, evs=regimen)
    return simobs(model, sub, mod, obstimes=0:0.1:12)
end

mapslices(Foo, lpars, [1000,2])

a = reshape(Vector(1:16),(8,2))


df_temp = @from i in obs begin
    @where i.time > 7
    @select {i.time, i.obs, i.sd}
    @collect DataFrame
end

@linq obs |>
@select(:obs) |>
@transform(Speed = ifelse.(:obs > 2, :obs-2, 0))
names(obs)
describe(obs[:obs])

using RCall
