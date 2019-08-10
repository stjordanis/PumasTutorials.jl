using Pumas, LinearAlgebra

model1 = @model begin
    @param begin
        θ = VectorDomain(3)
        Ω = PDiagDomain(3)
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @pre begin
        CL = θ[1] * exp(η[1])
        V = θ[2] * exp(η[2])
    end

    @dynamics ImmediateAbsorptionModel

    @derived begin
        cp = @. 1000*(Central / V)
    end

end

ev = DosageRegimen(2000, time=0, rate=-2)

pop =  Population(map(i -> Subject(id=i, evs=ev),1:10))
