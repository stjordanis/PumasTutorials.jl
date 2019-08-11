using Pumas, LinearAlgebra, Plots, Query

iv = @model begin
  @param   begin
    tvcl ∈ RealDomain(lower=0)
    tvv ∈ RealDomain(lower=0)
    Ω ∈ PDiagDomain(2)
    σ_prop ∈ RealDomain(lower=0)
  end

  @random begin
    η ~ MvNormal(Ω)
  end

  @pre begin
    CL = tvcl * (Wt/70)^0.75 * exp(η[1])
    V  = tvv * (Wt/70) * exp(η[2])
  end

  @covariates Wt

  @dynamics ImmediateAbsorptionModel

  @derived begin
      cp := @. (Central / V)
      dv ~ @. Normal(cp, sqrt(cp^2*σ_prop))
    end
end

ev = DosageRegimen(2000, time=0)
choose_covariates() = (wt = rand(55:80))
pop =  Population(map(i -> Subject(id=i, evs=ev, cvs = (Wt = choose_covariates(),)),1:24))
param = (
  tvcl = 11.5,
  tvv  = 50,
  Ω = Diagonal([0.04,0.04]),
  σ_prop = 0.001
  )

sd_obstimes = [0, 0.25, 0.5, 0.75, 1, 2, 4, 8,
                    12, 16, 20, 21, 24]
obs = simobs(iv, pop, param, obstimes=sd_obstimes)
iv_sim = DataFrame(obs)|>
          @filter(!(_.time == 0 && _.evid == 0)) |>
          @mutate(cmt = 1) |>
          @mutate(route = "iv") |>
          DataFrame

timeu = u"hr"
concu = u"mg/L"
amtu  = u"mg"

iv_nca = read_nca(iv_sim, id=:id, time=:time, conc=:dv, amt=:amt,
    route=:route,timeu=timeu, concu=concu, amtu=amtu,llq=0.4concu)

NCA.auc(iv_nca)
NCA.auc_extrap_percent(iv_nca)
report = NCAReport(iv_nca)

iv_sim_7percBLQ = iv_sim |>
                  @filter()
sd_obstimes_7percBLQ = [0, 0.25, 0.5, 0.75, 1, 2, 4, 8, 12, 16, 20]
sd_obstimes_18percBLQ = [0, 0.25, 0.5, 0.75, 1, 2, 4, 8, 12, 16, 20, 21, 24]

plot(obs,  yscale = :log10)
