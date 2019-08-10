## Enterohepatic Recycling Model

```julia
ehc = @model begin
  @param begin
    θ ∈ VectorDomain(8)
    Ω ∈ PDiagDomain(4)
  end

  @random begin
    η ~ MvNormal(Ω)
  end

  @covariates Wt

  @pre begin
    CL = θ[1]*(Wt/11)^0.75*exp(η[1])
    V  = θ[2]*(Wt/11)*exp(η[2])
    duration =  (Central = θ[3],)
  end

  @dynamics begin
    Central' =  - Central*CL/V
  end

  @derived begin
    conc = Central/V
  end

end
```

and below are the set of parameters for this model

```julia
param = (
  θ = [0.792, 13.7, 5],
  Ω =Diagonal([0.04,0.04,0.04]),
  )
```

We can now use the data and model with parameters to simulate a zero-order absorption profile

```julia
sd = DosageRegimen(100, time=0, rate=-2)
choose_covariates() = (Wt = rand(55:80),dose=100)
sd_with_covariates = Population(map(i -> Subject(id=i,evs=sd,cvs=choose_covariates()),1:10))
```

Note that while modeling duration, `rate` data item has to be set to `< 0`

```julia
obs = simobs(zoabs, sd_with_covariates, param, obstimes=0:0.1:48)
plot(obs, title =  "Zero-order absorption")
```

## Erlang absorption models with two compartment dsitribution

```julia
erlangabs = @model begin
  @param begin
    θ ∈ VectorDomain(5)
    Ω ∈ PDiagDomain(4)
  end

  @random begin
    η ~ MvNormal(Ω)
  end

  @covariates Wt dose

  @pre begin
    CL = θ[1]*(Wt/70)^0.75*exp(η[1])
    V1  = θ[2]*(Wt/70)*exp(η[2])
    Ktr = θ[3]*exp(η[3])
    Q  = θ[4]*(Wt/70)^0.75*exp(η[4])
    V2  = θ[5]*(Wt/70)
    K  = CL/V1
    K21 = Q/V2
    K12 = Q/V1
    N = 5 # number of compartments
  end

  @dynamics begin
    Central' = Ktr*dose*((Ktr*t)^N/factorial(N))*exp(-Ktr*t) - Central*K + K21*Periph -K12*Central
    Periph' = - K21*Periph + K12*Central
  end

  @derived begin
    conc = Central/V1
  end

end
```

```julia
sd = DosageRegimen(100, time=0, rate=0)
choose_covariates() = (Wt = rand(55:80),dose=100)
sd_with_covariates = Population(map(i -> Subject(id=i,evs=sd,cvs=choose_covariates()),1:10))
```

```julia
param = (
  θ = [26.3, 75.9, 4, 23.9, 43.1],
  Ω = Diagonal([0.04,0.04,0.04,0.04]),
  σ_prop = 0
  )
```
We can now use the data and model with parameters to simulate a Erlang absorption profile

```julia
obs = simobs(erlangabs, sd_with_covariates, param, obstimes=0:0.1:24)
plot(obs)
```


# 5 compartment model with Transit absorption for single dose

```julia
transabs = @model begin
  @param begin
    θ ∈ VectorDomain(6)
    Ω ∈ PDiagDomain(3)
  end

  @random begin
    η ~ MvNormal(Ω)
  end

  @covariates Wt dose

  @pre begin
    CL = θ[1]*(Wt/70)^0.75*exp(η[1])
    V  = θ[2]*(Wt/70)*exp(η[2])
    BIO = θ[3] # the dose was corrected by bioavailability already
    Ka = θ[4]
    MTT =θ[5]*exp(η[3])
    NN = θ[6]
    K  = CL/V
    Ktr = (NN+1)/MTT
    LNFAC=log(2.5066)+(NN+0.5)*log(NN)-NN
  end

  @dynamics begin
    Absorption' = exp(log(BIO*dose+.00001)+log(Ktr)+NN*log(Ktr*t+0.00001)-Ktr*t-LNFAC)-Ka*Absorption
    Central' = Ka*Absorption-Central*K
  end

  @derived begin
    conc = Central/V
  end

end
```

```julia
param = (
  θ = [3.54, 487, 1, 1.77, 6, 14],
  Ω = Diagonal([0.0361,0.16,0.27])
  )
```
We can now use the data and model with parameters to simulate a transit absorption profile

```julia
obs = simobs(transabs, sd_with_covariates, param, obstimes=0:0.1:120)
plot(obs)
```


$PROB Enterohepatic recycling
$INPUT ID TIME DV AMT CMT RATE EVID MDV
$DATA EHR_DATA.csv
IGNORE=#
$SUBROUTINE ADVAN6 TOL=4
$MODEL NCOMP=3
COMP=(CENTRAL,DEFOBS,DEFDOSE)
COMP=(PERIPH)
COMP=(ACCUM)
$PK
;duration (h) of the zero-order process
D1=THETA(1)*EXP(ETA(1))
;lag time for the zero-order process
ALAG1=THETA(2)
;disposition parameters
CL=THETA(3)*EXP(ETA(2))
V1=THETA(4)*EXP(ETA(3))
Q2=THETA(5)*EXP(ETA(4))
V2=THETA(6)
;Enterohepatic recycling parameters
K13=THETA(7)
T31=THETA(8)
S1=V1/1000
K10=CL/V1
K12=Q2/V1
K21=Q2/V2
;times for release from the accumulation compartment
MEA1=4+T31
MEA2=9+T31
MEA3=23+T31
$DES
; Default situation, only accumulation in the compt 3
DADT(1) = K21*A(2)-(K10+K12+K13)*A(1)
DADT(2) = K12*A(1)-K21*A(2)
DADT(3) = K13*A(1)
; First release
IF (TIME.GT.4.AND.TIME.LE.MEA1) THEN
DADT(1) = K21*A(2)-(K10+K12+K13)*A(1)+A(3)/T31
DADT(3) = K13*A(1)-A(3)/T31
ENDIF
; Second release
IF (TIME.GT.9.AND.TIME.LE.MEA2) THEN
DADT(1) = K21*A(2)-(K10+K12+K13)*A(1)+A(3)/T31
DADT(3) = K13*A(1)-A(3)/T31
ENDIF
; Third release
IF (TIME.GT.23.AND.TIME.LE.MEA3) THEN
DADT(1) = K21*A(2)-(K10+K12+K13)*A(1)+A(3)/T31
DADT(3) = K13*A(1)-A(3)/T31
ENDIF
$ERROR
A1=A(1)*1000
A2=A(2)*1000
A3=A(3)*1000
IPRED=F
;combined error
Y=F*(1+ERR(1))+ERR(2)
$THETA
0.196 ;1 D1
0.123 ;2 Alag1
0.518 ;3 CL
31.1 ;4 V1
2.23 ;5 Q2
47.1 ;6 V2
0.024 ;7 K13
1.55 ;8 T31
$SIGMA
0.0402
0.220
$OMEGA
0.500
0.312
0.103
0.584
ENTEROHEPATIC RECYCLING MODEL 381
382 POPULATION MODELS FOR DRUG ABSORPTION AND ENTEROHEPATIC RECYCLING
$SIMULATION (123456789) ONLYSIM
$TABLE NOPRINT ONEHEADER FILE=EHR.PAR
ID TIME DV AMT EVID MDV A1 A2 A3
ALAG1 D1 CL V1 V2 Q2 K13 T31 IPRED
