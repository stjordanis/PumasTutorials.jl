using PuMaS
using LinearAlgebra
using Plots


ir1 = @model begin
    @param begin
        θ ∈ VectorDomain(6)
        Ω ∈ PSDDomain(5)
    end

    @random begin
        η ~ MvNormal(Ω)
    end 


    @pre begin

        Ka1     = θ[1]
        CL      = θ[2]*exp(η[1])
        Vc      = θ[3]*exp(η[2])
        Kin     = θ[4]*exp(η[3])
        Kout    = θ[5]*exp(η[4])
        IC50    = θ[6]*exp(η[5])
        IMAX    = 1
    end


    @init begin
        Resp = Kin/Kout

    end


    @dynamics begin
        Ev1'    = -Ka1*Ev1
        Cent'   =  Ka1*Ev1 - (CL/Vc)*Cent
        Resp'   =  Kin*(1-(IMAX*(Cent/Vc)/(IC50+(Cent/Vc)))) - Kout*Resp
    end


    @derived begin
        cp     = Cent / Vc
        resp   = Resp
    end
end




regimen = DosageRegimen(1500, time=0,cmt=1)


subject = Subject(id=1,evs=regimen)


p = (θ = [
          0.5, # Ka1  Absorption rate constant 1 (1/time)
          5, # CL   Clearance (volume/time)
          30, # Vc   Central volume (volume)
          20, # Kin  Response in rate constant (1/time)
          0.01, # Kout Response out rate constant (1/time)
          10, # IC50 Concentration for 50% of max inhibition (mass/volume)
          1, # IMAX Maximum inhibition
          ]),
       (Ω = PDMat(diagm(0 => [0.01,0.01,0.01,0.01,0.01]))
      )


sim = simobs(ir1, subject, p, obstimes=0:0.1:200, parallel_type=PuMaS.Serial, alg=Rodas5())


plot(sim)
