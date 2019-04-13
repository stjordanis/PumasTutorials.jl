using PuMaS, LinearAlgebra, Plots
using CSV, DataFrames, Query
using GLM, PuMaS.NCA, DiffEqSensitivity, Random
using DataStructures

model = @model begin

    @param begin
        θ ∈ VectorDomain(11)
    end

    @pre begin

        fup = 0.42
        fumic = 0.711
        WEIGHT = 73
        MPPGL = 30.3
        MPPGI = 0
        C_OUTPUT = 6.5
        VmaxH = 40
        VmaxG = 40
        KmH = 9.3
        KmG = 9.3
        BP = 1

        Kpad = θ[1]
        Kpbo = θ[2]
        Kpbr = θ[3]
        Kpgu = θ[4]
        Kphe = θ[5]
        Kpki = θ[6]
        Kpli = θ[7]
        Kplu = θ[8]
        Kpmu = θ[9]
        Kpsp = θ[10]
        Kpre = θ[11]

        MW = 349.317
        logP = 2.56
        S_lumen = 0.39*1000
        L = 280
        d = 2.5
        PF = 1.57
        VF = 6.5
        MF = 13
        ITT = 3.32
        A = 7440
        B = 1e7
        alpha = 0.6
        beta = 4.395
        fabs = 1
        fdis = 1
        fperm = 1

        Vad = 18.2
        Vbo = 10.5
        Vbr = 1.45
        VguWall = 0.65
        VguLumen = 0.35
        Vgu = VguWall + VguLumen
        Vhe = 0.33
        Vki = 0.31
        Vli = 1.8
        Vlu = 0.5
        Vmu = 29
        Vsp = 0.15
        Vbl = 5.6
        Vve = 0.705*Vbl
        Var = 0.295*Vbl
        Vre = WEIGHT - (Vli+Vki+Vsp+Vhe+Vlu+Vbo+Vbr+Vmu+Vad+VguWall+Vbl)

        FQad = 0.05
        FQbo = 0.05
        FQbr = 0.12
        FQgu = 0.16
        FQhe = 0.04
        FQki = 0.19
        FQli = 0.255
        FQmu = 0.17
        FQsp = 0.03

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
        GUTLUMEN' = -kd*VguLumen*(f*(GUTLUMEN/VguLumen) + (1-f)*S_lumen) - kt*GUTLUMEN
        GUTWALL' = kd*VguLumen*(f*(GUTLUMEN/VguLumen) + (1-f)*S_lumen) - ka*GUTWALL - CLintGut*(GUTWALL/VguWall)
        GUT' = ka*GUTWALL + Qgu*((ART/Var) - (GUT/VguWall)/(Kpgu/BP))
        ADIPOSE' = Qad*((ART/Var) - (ADIPOSE/Vad)/(Kpad/BP))
        BRAIN' = Qbr*((ART/Var) - (BRAIN/Vbr)/(Kpbr/BP))
        HEART' = Qhe*((ART/Var) - (HEART/Vhe)/(Kphe/BP))
        KIDNEY' = Qki*((ART/Var) - (KIDNEY/Vki)/(Kpki/BP)) - CLrenal*(((KIDNEY/Vki)*fup)/(Kpki/BP))
        LIVER' = Qgu*((GUT/VguWall)/(Kpgu/BP)) + Qsp*((SPLEEN/Vsp)/(Kpsp/BP)) + Qha*(ART/Var) - Qli*((LIVER/Vli)/(Kpli/BP)) - CLintHep*(((LIVER/Vli)*fup)/(Kpli/BP))
        LUNG' = Qlu*((VEN/Vve) - (LUNG/Vlu)/(Kplu/BP))
        MUSCLE' = Qmu*((ART/Var) - (MUSCLE/Vmu)/(Kpmu/BP))
        SPLEEN' = Qsp*((ART/Var) - (SPLEEN/Vsp)/(Kpsp/BP))
        BONE' = Qbo*((ART/Var) - (BONE/Vbo)/(Kpbo/BP))
        REST' = Qre*((ART/Var) - (REST/Vre)/(Kpre/BP))
        VEN' = Qad*((ADIPOSE/Vad)/(Kpad/BP)) + Qbr*((BRAIN/Vbr)/(Kpbr/BP)) + Qhe*((HEART/Vhe)/(Kphe/BP)) + Qki*((KIDNEY/Vki)/(Kpki/BP)) + Qli*((LIVER/Vli)/(Kpli/BP)) + Qmu*((MUSCLE/Vmu)/(Kpmu/BP)) + Qbo*((BONE/Vbo)/(Kpbo/BP)) + Qre*((REST/Vre)/(Kpre/BP)) - Qlu*(VEN/Vve)
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
    end
end

type = 3
logP = 2.56
pKa = 1.76
fup = .42
BP = 1
BW = 73


regimen = DosageRegimen(292, time = 0, addl=13, ii=12, cmt=14, rate = 292, ss = 1)

sub1 = Subject(id=1,evs=regimen)

p = (θ = [
9.89, #Kpad
7.91, #Kpbo
7.35, #Kpbr
5.82, #Kpgu
1.95, #Kphe
2.9, #Kpki
4.66, #Kpli
0.83, #Kplu
2.94, #Kpmu
2.96, #Kpsp
4 #Kpre
],)

pred_PT = simobs(model, sub1, p, obstimes=0:0.1:12)


p_Berez = (θ = [
10.7, #Kpad
7.33, #Kpbo
6.87, #Kpbr
5.82, #Kpgu ??
1.92, #Kphe
2.8, #Kpki
4.4, #Kpli
0.91, #Kplu
2.83, #Kpmu
2.85, #Kpsp
4 #Kpre
],)

pred_Berez = simobs(model, sub1, p_Berez, obstimes=0:0.1:12) #slightly different numbers - could be due to Kpgu?

p_RR = (θ = [
11.5, #Kpad
2.86, #Kpbo
6.34, #Kpbr
6.69, #Kpgu
2.97, #Kphe
3.26, #Kpki
3.55, #Kpli
4.27, #Kplu
2.17, #Kpmu
2.02, #Kpsp
4 #Kpre
],)

pred_RR = simobs(model, sub1, p_RR, obstimes = 0:0.1:12) #slightly different numbers

p_Schmitt = (θ = [
141, #Kpad
3.88, #Kpbo
17.4, #Kpbr
5.82, #Kpgu
17.7, #Kphe
10.2, #Kpki
13.1, #Kpli
6.79, #Kplu
2.5, #Kpmu
4.12, #Kpsp
4 #Kpre
],)

pred_Schmitt = simobs(model, sub1, p_Schmitt, obstimes = 0:0.1:12) #slightly different numbers


p_pksim = (θ = [
122, #Kpad
42, #Kpbo
17.4, #Kpbr
5.82, #Kpgu
16.1, #Kphe
8.84, #Kpki
11.5, #Kpli
2.49, #Kplu
2.93, #Kpmu
3.43, #Kpsp
4 #Kpre
],)

pred_pksim = simobs(model, sub1, p_pksim, obstimes = 0:0.1:12) #slightly different numbers


#csv import attempt
obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/obs.csv"))
ZT = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/ZT.csv"))


#plot
plot(pred_PT.times, pred_PT[:Cvenn], linecolor = :black, label = "PT", title = "Adult 4 mg/kg IV")
plot!(pred_Berez.times, pred_Berez[:Cvenn], line = [:dash], linecolor = :black, label = "Berez")
plot!(pred_RR.times, pred_RR[:Cvenn], line = [:dot], linecolor = :black, label = "RR")
plot!(pred_Schmitt.times, pred_Schmitt[:Cvenn], line = :dashdot, linecolor = :black, label = "Schmitt")
plot!(pred_pksim.times, pred_pksim[:Cvenn], line = :dashdotdot, linecolor = :black, label = "PK Sim")
plot!(obs.time, obs[:obs], seriestype=:scatter, yerror = obs[:sd], markercolor = :black, seriescolor = :black, label = "observed")
plot!(ZT.time, ZT[:ZT], linecolor = :lightgray, label = "ZT")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")


#Table 1: This chunk tables out the prediction errors of different calculation methods
#pad observed data till 12 hours
df_temp = @from i in obs begin
    @where i.time > 7
    @select {i.time, i.obs, i.sd}
    @collect DataFrame
end

mod = coef(lm(@formula(obs ~ time), df_temp))

st = DataFrame(time = [12])
st2 = predict((lm(@formula(obs ~ time), df_temp)), st)

df_temp = DataFrame(Column1 = 14, time = 12, obs = st2, sd = missing)
append!(obs, df_temp)

## calculate AUCs for all methods
auc_obs = NCA.auc(obs[:obs], obs[:time], interval = (0, 12))
auc_PT = NCA.auc(pred_PT[:Cvenn], pred_PT.times, interval = (0,12))
auc_Berez = NCA.auc(pred_Berez[:Cvenn], pred_Berez.times, interval = (0,12))
auc_RR = NCA.auc(pred_RR[:Cvenn], pred_RR.times, interval = (0,12))
auc_Schmitt = NCA.auc(pred_Schmitt[:Cvenn], pred_Schmitt.times, interval = (0,12))
auc_pksim = NCA.auc(pred_pksim[:Cvenn], pred_pksim.times, interval = (0,12))

cmax_obs = NCA.cmax(obs[:obs], obs[:time])
cmax_PT = NCA.cmax(pred_PT[:Cvenn], pred_PT.times)
cmax_Berez = NCA.cmax(pred_Berez[:Cvenn], pred_Berez.times)
cmax_RR = NCA.cmax(pred_RR[:Cvenn], pred_RR.times)
cmax_Schmitt = NCA.cmax(pred_Schmitt[:Cvenn], pred_Schmitt.times)
cmax_pksim = NCA.cmax(pred_pksim[:Cvenn], pred_pksim.times)

#need to find out how to mutate


## Figure 4a; Model 1 with 4 mg/kg IV infusion
regimen2 = DosageRegimen(292, time = 0, addl=7, ii=12, cmt=14, rate = 292, ss = 1)

sub2 = Subject(id=1,evs=regimen2)
sim = simobs(model, sub2, p_RR, obstimes = 0:0.1:12) #slightly different numbers

#csv import attempt
obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/obs.csv"))
ZT = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/ZT.csv"))

plot(sim.times, sim[:Cvenn], linecolor = :black, label = "model 1", title = "a Adult 4 mg/kg IV")
plot!(obs.time, obs[:obs], seriestype=:scatter, yerror = obs[:sd], markercolor = :black, seriescolor = :black, label = "observed")
plot!(ZT.time, ZT[:ZT], linecolor = :lightgray, label = "ZT")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")


## Figure 3b; Model 2 with 4 mg/kg IV infusion
obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4a_obs.csv"))
ZT = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4a_ZT.csv"))


model2 = @model begin
    @param begin
        fup ∈ ConstDomain(0.42)
        fumic ∈ ConstDomain(0.711)
        WEIGHT ∈ ConstDomain(19)
        MPPGL ∈ ConstDomain(26)
        MPPGI ∈ ConstDomain(0)
        C_OUTPUT ∈ ConstDomain(3.4)
        VmaxH ∈ ConstDomain(120.5)
        VmaxG ∈ ConstDomain(120.5)
        KmH ∈ ConstDomain(11)
        KmG ∈ ConstDomain(11)
        BP ∈ ConstDomain(1)

        Kpad ∈ ConstDomain(9.89)
        Kpbo ∈ ConstDomain(7.91)
        Kpbr ∈ ConstDomain(7.35)
        Kpgu ∈ ConstDomain(5.82)
        Kphe ∈ ConstDomain(1.95)
        Kpki ∈ ConstDomain(2.9)
        Kpli ∈ ConstDomain(4.66)
        Kplu ∈ ConstDomain(0.83)
        Kpmu ∈ ConstDomain(2.94)
        Kpsp ∈ ConstDomain(2.96)
        Kpre ∈ ConstDomain(4)

        MW ∈ ConstDomain(349.317)
        logP ∈ ConstDomain(2.56)
        S_lumen ∈ ConstDomain(0.39*1000)
        L ∈ ConstDomain(170)
        d ∈ ConstDomain(2.5)
        PF ∈ ConstDomain(1.57)
        VF ∈ ConstDomain(6.5)
        MF ∈ ConstDomain(13)
        ITT ∈ ConstDomain(3.32)
        A ∈ ConstDomain(7440)
        B ∈ ConstDomain(1e7)
        alpha ∈ ConstDomain(0.6)
        beta ∈ ConstDomain(4.395)
        fabs ∈ ConstDomain(1)
        fdis ∈ ConstDomain(1)
        fperm ∈ ConstDomain(1)

        Vad ∈ ConstDomain(5.5)
        Vbo ∈ ConstDomain(2.43)
        Vbr ∈ ConstDomain(1.31)
        VguWall ∈ ConstDomain(0.22)
        VguLumen ∈ ConstDomain(0.117)

        Vhe ∈ ConstDomain(0.085)
        Vki ∈ ConstDomain(0.11)
        Vli ∈ ConstDomain(0.467)
        Vlu ∈ ConstDomain(0.125)
        Vmu ∈ ConstDomain(5.6)
        Vsp ∈ ConstDomain(0.05)
        Vbl ∈ ConstDomain(1.5)

        FQad ∈ ConstDomain(0.05)
        FQbo ∈ ConstDomain(0.05)
        FQbr ∈ ConstDomain(0.12)
        FQgu ∈ ConstDomain(0.16)
        FQhe ∈ ConstDomain(0.04)
        FQki ∈ ConstDomain(0.19)
        FQli ∈ ConstDomain(0.255)
        FQmu ∈ ConstDomain(0.17)
        FQsp ∈ ConstDomain(0.03)
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
        GUTLUMEN' = -kd*VguLumen*(f*(GUTLUMEN/VguLumen) + (1-f)*S_lumen) - kt*GUTLUMEN
        GUTWALL' = kd*VguLumen*(f*(GUTLUMEN/VguLumen) + (1-f)*S_lumen) - ka*GUTWALL - CLintGut*(GUTWALL/VguWall)
        GUT' = ka*GUTWALL + Qgu*((ART/Var) - (GUT/VguWall)/(Kpgu/BP))
        ADIPOSE' = Qad*((ART/Var) - (ADIPOSE/Vad)/(Kpad/BP))
        BRAIN' = Qbr*((ART/Var) - (BRAIN/Vbr)/(Kpbr/BP))
        HEART' = Qhe*((ART/Var) - (HEART/Vhe)/(Kphe/BP))
        KIDNEY' = Qki*((ART/Var) - (KIDNEY/Vki)/(Kpki/BP)) - CLrenal*(((KIDNEY/Vki)*fup)/(Kpki/BP))
        LIVER' = Qgu*((GUT/VguWall)/(Kpgu/BP)) + Qsp*((SPLEEN/Vsp)/(Kpsp/BP)) + Qha*(ART/Var) - Qli*((LIVER/Vli)/(Kpli/BP)) - CLintHep*(((LIVER/Vli)*fup)/(Kpli/BP))
        LUNG' = Qlu*((VEN/Vve) - (LUNG/Vlu)/(Kplu/BP))
        MUSCLE' = Qmu*((ART/Var) - (MUSCLE/Vmu)/(Kpmu/BP))
        SPLEEN' = Qsp*((ART/Var) - (SPLEEN/Vsp)/(Kpsp/BP))
        BONE' = Qbo*((ART/Var) - (BONE/Vbo)/(Kpbo/BP))
        REST' = Qre*((ART/Var) - (REST/Vre)/(Kpre/BP))
        VEN' = Qad*((ADIPOSE/Vad)/(Kpad/BP)) + Qbr*((BRAIN/Vbr)/(Kpbr/BP)) + Qhe*((HEART/Vhe)/(Kphe/BP)) + Qki*((KIDNEY/Vki)/(Kpki/BP)) + Qli*((LIVER/Vli)/(Kpli/BP)) + Qmu*((MUSCLE/Vmu)/(Kpmu/BP)) + Qbo*((BONE/Vbo)/(Kpbo/BP)) + Qre*((REST/Vre)/(Kpre/BP)) - Qlu*(VEN/Vve)
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
    end
end

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

p_RR_ped = (fixeffs_ped..., Kpad = 11.5, Kpbo = 2.86, Kpbr = 6.34, Kpgu = 6.69, Kphe = 2.97,
        Kpki = 3.26, Kpli = 3.55, Kplu = 4.27, Kpmu = 2.17, Kpsp = 2.02, Kpre = 4)

regimen3 = DosageRegimen(76, time = 0, addl=7, ii=12, cmt=14, rate = 57, ss = 1)
sub3 = Subject(id=1,evs=regimen3)
sim = simobs(model2, sub3, p_RR_ped, obstimes = 0:0.1:12)


plot(sim.times, sim[:Cvenn], linecolor = :black, label = "model 2", title = "b Pediatric 4 mg/kg IV")
plot!(obs.time, obs[:obs], seriestype=:scatter, yerror = obs[:sd], markercolor = :black, seriescolor = :black, label = "observed")
plot!(ZT.time, ZT[:ZT], linecolor = :lightgray, label = "ZT")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")

#Figure 3c; Model 1 with 200 mg PO
obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig3b_obs.csv"))
ZT = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig3b_ZT.csv"))

regimen4 = DosageRegimen(200, time = 0, addl=7, ii=12, cmt=1, ss = 1)
sub4 = Subject(id=1, evs=regimen4)
sim = simobs(model, sub4, p_RR, obstimes = 0:0.1:12)

plot(sim.times, sim[:Cvenn], linecolor = :black, label = "model 1", title = "c Adult 200 mg PO")
plot!(obs.time, obs[:obs], seriestype=:scatter, yerror = obs[:sd], markercolor = :black, seriescolor = :black, label = "observed")
plot!(ZT.time, ZT[:ZT], linecolor = :lightgray, label = "ZT")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")


#Figure 3d; Model 2 with 4 mg/kg PO
obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4d_obs.csv"))
ZT = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4d_ZT.csv"))

regimen5 = DosageRegimen(76, time = 0, addl=7, ii=12, cmt=1, ss = 1)
sub5 = Subject(id=1, evs=regimen5)
sim = simobs(model2, sub5, p_RR_ped, obstimes = 0:0.1:12)

plot(sim.times, sim[:Cvenn], linecolor = :black, label = "model 1", title = "d Pediatric 4 mg/kg PO")
plot!(obs.time, obs[:obs], seriestype=:scatter, yerror = obs[:sd], markercolor = :black, seriescolor = :black, label = "observed")
plot!(ZT.time, ZT[:ZT], linecolor = :lightgray, label = "ZT")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")


### sensitivity analyses for absorption parameters: permeability, intestinal transit time and solubility ###

model = @model begin

    @param begin
        fup ∈ ConstDomain(0.42)
        fumic ∈ ConstDomain(0.711)
        WEIGHT ∈ ConstDomain(73)
        MPPGL ∈ ConstDomain(30.3)
        MPPGI ∈ ConstDomain(0)
        C_OUTPUT ∈ ConstDomain(6.5)
        VmaxH ∈ ConstDomain(40)
        VmaxG ∈ ConstDomain(40)
        KmH ∈ ConstDomain(9.3)
        KmG ∈ ConstDomain(9.3)
        BP ∈ ConstDomain(1)

        Kpad ∈ ConstDomain(9.89)
        Kpbo ∈ ConstDomain(7.91)
        Kpbr ∈ ConstDomain(7.35)
        Kpgu ∈ ConstDomain(5.82)
        Kphe ∈ ConstDomain(1.95)
        Kpki ∈ ConstDomain(2.9)
        Kpli ∈ ConstDomain(4.66)
        Kplu ∈ ConstDomain(0.83)
        Kpmu ∈ ConstDomain(2.94)
        Kpsp ∈ ConstDomain(2.96)
        Kpre ∈ ConstDomain(4)


        MW ∈ ConstDomain(349.317)
        logP ∈ ConstDomain(2.56)
        S_lumen ∈ ConstDomain(0.39*1000)
        L ∈ ConstDomain(280)
        d ∈ ConstDomain(2.5)
        PF ∈ ConstDomain(1.57)
        VF ∈ ConstDomain(6.5)
        MF ∈ ConstDomain(13)
        ITT ∈ ConstDomain(3.32)
        A ∈ ConstDomain(7440)
        B ∈ ConstDomain(1e7)
        alpha ∈ ConstDomain(0.6)
        beta ∈ ConstDomain(4.395)
        fabs ∈ ConstDomain(1)
        fdis ∈ ConstDomain(1)
        fperm ∈ ConstDomain(1)

        Vad ∈ ConstDomain(18.2)
        Vbo ∈ ConstDomain(10.5)
        Vbr ∈ ConstDomain(1.45)
        VguWall ∈ ConstDomain(0.65)
        VguLumen ∈ ConstDomain(0.35)

        Vhe ∈ ConstDomain(0.33)
        Vki ∈ ConstDomain(0.31)
        Vli ∈ ConstDomain(1.8)
        Vlu ∈ ConstDomain(0.5)
        Vmu ∈ ConstDomain(29)
        Vsp ∈ ConstDomain(0.15)
        Vbl ∈ ConstDomain(5.6)

        FQad ∈ ConstDomain(0.05)
        FQbo ∈ ConstDomain(0.05)
        FQbr ∈ ConstDomain(0.12)
        FQgu ∈ ConstDomain(0.16)
        FQhe ∈ ConstDomain(0.04)
        FQki ∈ ConstDomain(0.19)
        FQli ∈ ConstDomain(0.255)
        FQmu ∈ ConstDomain(0.17)
        FQsp ∈ ConstDomain(0.03)
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
        GUTLUMEN' = -kd*VguLumen*(f*(GUTLUMEN/VguLumen) + (1-f)*S_lumen) - kt*GUTLUMEN
        GUTWALL' = kd*VguLumen*(f*(GUTLUMEN/VguLumen) + (1-f)*S_lumen) - ka*GUTWALL - CLintGut*(GUTWALL/VguWall)
        GUT' = ka*GUTWALL + Qgu*((ART/Var) - (GUT/VguWall)/(Kpgu/BP))
        ADIPOSE' = Qad*((ART/Var) - (ADIPOSE/Vad)/(Kpad/BP))
        BRAIN' = Qbr*((ART/Var) - (BRAIN/Vbr)/(Kpbr/BP))
        HEART' = Qhe*((ART/Var) - (HEART/Vhe)/(Kphe/BP))
        KIDNEY' = Qki*((ART/Var) - (KIDNEY/Vki)/(Kpki/BP)) - CLrenal*(((KIDNEY/Vki)*fup)/(Kpki/BP))
        LIVER' = Qgu*((GUT/VguWall)/(Kpgu/BP)) + Qsp*((SPLEEN/Vsp)/(Kpsp/BP)) + Qha*(ART/Var) - Qli*((LIVER/Vli)/(Kpli/BP)) - CLintHep*(((LIVER/Vli)*fup)/(Kpli/BP))
        LUNG' = Qlu*((VEN/Vve) - (LUNG/Vlu)/(Kplu/BP))
        MUSCLE' = Qmu*((ART/Var) - (MUSCLE/Vmu)/(Kpmu/BP))
        SPLEEN' = Qsp*((ART/Var) - (SPLEEN/Vsp)/(Kpsp/BP))
        BONE' = Qbo*((ART/Var) - (BONE/Vbo)/(Kpbo/BP))
        REST' = Qre*((ART/Var) - (REST/Vre)/(Kpre/BP))
        VEN' = Qad*((ADIPOSE/Vad)/(Kpad/BP)) + Qbr*((BRAIN/Vbr)/(Kpbr/BP)) + Qhe*((HEART/Vhe)/(Kphe/BP)) + Qki*((KIDNEY/Vki)/(Kpki/BP)) + Qli*((LIVER/Vli)/(Kpli/BP)) + Qmu*((MUSCLE/Vmu)/(Kpmu/BP)) + Qbo*((BONE/Vbo)/(Kpbo/BP)) + Qre*((REST/Vre)/(Kpre/BP)) - Qlu*(VEN/Vve)
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
    end
end

fixeffs = (fup = 0.42, fumic = 0.711, WEIGHT = 73, MPPGL = 30.3, MPPGI = 0,
    C_OUTPUT = 6.5, VmaxH = 40, VmaxG = 40, KmH = 9.3, KmG = 9.3, BP = 1,
    Kpad = 9.89, Kpbo = 7.91, Kpbr = 7.35, Kpgu = 5.82, Kphe = 1.95, Kpki = 2.9,
    Kpli = 4.66, Kplu = 0.83, Kpmu = 2.94, Kpsp = 2.96, Kpre = 4, MW = 349.317,
    logP = 2.56, S_lumen = 0.39*1000, L = 280, d = 2.5, PF = 1.57, VF = 6.5,
    MF = 13, ITT = 3.32, A = 7440, B = 1e7, alpha = 0.6, beta = 4.395, fabs = 1,
    fdis = 1, fperm = 1, Vad = 18.2, Vbo = 10.5, Vbr = 1.45, VguWall = 0.65,
    VguLumen = 0.35, Vhe = 0.33, Vki = 0.31, Vli = 1.8, Vlu = 0.5, Vmu = 29,
    Vsp = 0.15, Vbl = 5.6, FQad = 0.05, FQbo = 0.05, FQbr = 0.12, FQgu = 0.16,
    FQhe = 0.04, FQki = 0.19, FQli = 0.255, FQmu = 0.17, FQsp = 0.03)

regimen_s = DosageRegimen(200, time = 0, addl=13, ii=12, cmt=1, ss = 1)

sub_s = Subject(id=1,evs=regimen_s)

#changing S_lumen
_fixeffs = (fixeffs...,S_lumen = 0.39*500)
_fixeffs2 = (fixeffs...,S_lumen = 0.39*1000)
_fixeffs3 = (fixeffs...,S_lumen = 0.39*2000)

#Kpmu sensitivity analysis for abstract
_Kpmu1 = (fixeffs..., Kpmu = 1.5)
_Kpmu2 = (fixeffs..., Kpmu = 3)
_Kpmu3 = (fixeffs..., Kpmu = 6)

pred_kpmu1 = simobs(model, sub_s, _Kpmu1, obstimes=0:0.1:12)
pred_kpmu2 = simobs(model, sub_s, _Kpmu2, obstimes=0:0.1:12)
pred_kpmu3 = simobs(model, sub_s, _Kpmu3, obstimes=0:0.1:12)

p1 = plot(pred_kpmu1.times, pred_kpmu1[:Cvenn], linecolor = :black, label = "1.5", legendtitle = "Kpmu", legendfontsize = 12)
plot!(pred_kpmu2.times, pred_kpmu2[:Cvenn], linecolor = :black, line = :dash, label = "3")
plot!(pred_kpmu3.times, pred_kpmu3[:Cvenn], linecolor = :black, line = :dot, label = "6")
xlabel!("Time (h)")
ylabel!("Plasma concentration (mg/L)")
p1


#BP sensitivity analysis for abstract
_BP1 = (fixeffs..., BP = 0.5)
_BP2 = (fixeffs..., BP = 1)
_BP3 = (fixeffs..., BP = 2)

pred_BP1 = simobs(model, sub_s, _BP1, obstimes=0:0.1:12)
pred_BP2 = simobs(model, sub_s, _BP2, obstimes=0:0.1:12)
pred_BP3 = simobs(model, sub_s, _BP3, obstimes=0:0.1:12)

p2 = plot(pred_BP1.times, pred_BP1[:Cvenn], linecolor = :black, label = "0.5", legendtitle = "BP", legendfontsize = 12)
plot!(pred_BP2.times, pred_BP2[:Cvenn], linecolor = :black, line = :dash, label = "1")
plot!(pred_BP3.times, pred_BP3[:Cvenn], linecolor = :black, line = :dot, label = "2")
xlabel!("Time (h)")
ylabel!("Plasma concentration (mg/L)")


pred_1 = simobs(model, sub_s, _fixeffs, obstimes=0:0.1:12)
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

pred_4 = simobs(model, sub_s, _fixeffs4, obstimes=0:0.1:12)
pred_5 = simobs(model, sub_s, _fixeffs5, obstimes=0:0.1:12)
pred_6 = simobs(model, sub_s, _fixeffs6, obstimes=0:0.1:12)

plot(pred_4.times, pred_4[:Cvenn], linecolor = :black, label = "0.073", title = "Peff x 10-4 (cm/s)")
plot!(pred_5.times, pred_5[:Cvenn], linecolor = :black, line = :dash, label = "0.145")
plot!(pred_6.times, pred_6[:Cvenn], linecolor = :black, line = :dot, label = "0.29")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")

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



###Figure 7
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

obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig3b_obs.csv"))
ZT = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig3b_ZT.csv"))

regimen = DosageRegimen(200, time = 0, addl=7, ii=12, cmt=1, ss = 1)

sub = Subject(id=1,evs=regimen)

sim = simobs(model, sub, p_RR, obstimes=0:0.1:12)
sim_mppgi = simobs(model, sub, model3, obstimes=0:0.1:12)
sim_fperm = simobs(model, sub, model5, obstimes=0:0.1:12)
sim_mppgi_fperm = simobs(model, sub, model7, obstimes=0:0.1:12)

plot(sim.times, sim[:Cvenn], line = [:dash],linecolor = :black, label = "model 1", title = "Adult 200mg PO")
plot!(sim_mppgi.times, sim_mppgi[:Cvenn], line = [:dashdot], linecolor = :black, label = "model 3")
plot!(sim_fperm.times, sim_fperm[:Cvenn], line = [:dot], linecolor = :black, label = "model 5")
plot!(sim_mppgi_fperm.times, sim_mppgi_fperm[:Cvenn], linecolor = :black, label = "model 7")
plot!(obs.time, obs[:obs], seriestype=:scatter, yerror = obs[:sd], markercolor = :black, seriescolor = :black, label = "observed")
plot!(ZT.time, ZT[:ZT], linecolor = :lightgray, label = "ZT")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")


#Figure 5b; Model 4 with 4 mg/kg IV infusion
obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4d_obs.csv"))
ZT_Gu = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4e_ZT.csv"))
ZT = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4d_ZT.csv"))

regimen = DosageRegimen(76, time = 0, addl=7, ii=12, cmt=1, ss = 1)
sub = Subject(id=1,evs=regimen)

sim = simobs(model2, sub, p_RR_ped, obstimes=0:0.1:12)
sim_mppgi = simobs(model2, sub, model4, obstimes=0:0.1:12)
sim_fperm = simobs(model2, sub, model6, obstimes=0:0.1:12)
sim_mppgi_fperm = simobs(model2, sub, model8, obstimes=0:0.1:12)


plot(sim.times, sim[:Cvenn], line = [:dash],linecolor = :black, label = "model 1", title = "Pediatric 4mg/kg PO")
plot!(sim_mppgi.times, sim_mppgi[:Cvenn], line = [:dashdot], linecolor = :black, label = "model 4")
plot!(sim_fperm.times, sim_fperm[:Cvenn], line = [:dot], linecolor = :black, label = "model 6")
plot!(sim_mppgi_fperm.times, sim_mppgi_fperm[:Cvenn], linecolor = :black, label = "model 8")
plot!(obs.time, obs[:obs], seriestype=:scatter, yerror = obs[:sd], markercolor = :black, seriescolor = :black, label = "observed")
plot!(ZT.time, ZT[:ZT], linecolor = :lightgray, label = "ZT")
plot!(ZT_Gu.time, ZT_Gu[:ZT], linecolor = :lightgray, line = :dash, label ="ZT_Gu")
xlabel!("time (h)")
ylabel!("Plasma concentration (mg/L)")

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

df_temp = DataFrame(Column1 = 14, time = 12, obs = st2, sd = missing)
append!(obs, df_temp)


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

obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4a_obs.csv"))
ZT = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4a_ZT.csv"))

regimen = DosageRegimen(76, time = 0, addl=13, ii=12, cmt = 14, rate = 76, ss = 1)
sub = Subject(id=1, evs=regimen)

pred = simobs(model2, sub, p_RR_ped, obstimes=0:0.1:12)

auc_obs = NCA.auc(obs[:obs], obs[:time], interval = (0, 12))
auc_ZT = NCA.auc(ZT[:ZT], ZT[:time], interval = (0, 12))
auc_pred = NCA.auc(pred[:Cvenn], pred.times, interval = (0,12))

cmax_obs = NCA.cmax(obs[:obs], obs[:time])
cmax_ZT = NCA.cmax(ZT[:ZT], ZT[:time])
cmax_pred = NCA.cmax(pred[:Cvenn], pred.times)

#need to make into dataframe and then use mutate

#Table 3
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

df_temp = DataFrame(Column1 = 14, time = 12, obs = st2, sd = missing)
append!(obs, df_temp)

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

obs = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4d_obs.csv"))
ZT = DataFrame(CSV.File("/Users/praneeth/Documents/juliatrial/Fig4e_ZT.csv"))

regimen = DosageRegimen(76, time = 0, addl=13, ii=12, cmt = 1, ss = 1)
sub = Subject(id=1, evs=regimen)

pred_init = simobs(model2, sub, p_RR_ped, obstimes=0:0.1:12)
pred_mppgi = simobs(model2, sub, model4, obstimes=0:0.1:12)
pred_fperm = simobs(model2, sub, model6, obstimes=0:0.1:12)
pred_mppgi_fperm = simobs(model2, sub, model8, obstimes=0:0.1:12)

auc_obs = NCA.auc(obs[:obs], obs[:time], interval = (0, 12)) #auc value different by 0.7
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

#Figure 8
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
