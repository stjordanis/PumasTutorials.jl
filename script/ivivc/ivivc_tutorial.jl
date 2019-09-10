
using Pumas.IVIVC


vivo_data = read_vivo("./vivo_data.csv");


vivo_data = read_vivo("./vivo_data.csv", conc=:conc);


using CSV
df = CSV.read("./vivo_data.csv")
vivo_data = read_vivo(df);


vivo_fast_data = vivo_data[1]["fast"]


vivo_fast_data.time


get_avail_models()


vitro_data = read_vitro("./vitro_data.csv");
vitro_fast_data = vitro_data[1]["fast"]
estimate_fdiss(vitro_fast_data, :emax);


using Plots
scatter(vitro_fast_data.time, vitro_fast_data.conc, label="Observed data", xlabel="time", ylabel="Fdiss")
plot!(vitro_fast_data, label="Fitted model", legend=:topleft)


Emax_model(t, p) = @. (p[1] * (t ^ p[2])) / (p[3]^p[2] + t^p[2])
p0 = [vitro_fast_data.conc[end], 1.2, vitro_fast_data.time[2]]
lb = [0.0, 1.0, 0.0]
ub = [1.25, Inf, vitro_fast_data.time[end]]
estimate_fdiss(vitro_fast_data, Emax_model, p0=p0, lower_bound=lb, upper_bound=ub);


vitro_fast_data.pmin


estimate_fdiss(vitro_data[1]["medium"], :emax)
estimate_fdiss(vitro_data[1]["slow"], :emax);


uir_data = read_uir("./uir_data.csv")
estimate_uir(uir_data, :bateman);


estimate_uir(uir_data, :bateman, frac=1.0);


ka, kel, V = uir_data.pmin
@show ka
@show kel
@show V;


time = vivo_fast_data.time
conc = vivo_fast_data.conc
fabs, conc_auc = estimate_fabs(conc, time, kel, :wn)
@show fabs
@show conc_auc;


all_fabs_dict = Dict(); all_conc_auc = Dict()
for (formulation, profile) in vivo_data[1]
  all_fabs_dict[formulation], all_conc_auc[formulation] = estimate_fabs(profile.conc,
                                                              profile.time, kel, :wn)
end


ivivc_model = (form, time, x) -> x[1] * vitro_data[1][form](time * x[2])
# initial estimates of parameters, upper_bounds and lower_bounds
p0 = [0.8, 0.5]
ub = [1.25, 1.25]
lb = [0.0, 0.0];


mse(x, y) = sum(abs2.(x .- y))/length(x)
function errfun(x)
  err = 0.0
  for (form, prof) in vivo_data[1]   # (formulation => profile)
    err = err + mse(ivivc_model(form, prof.time, x), all_fabs_dict[form])
  end
  return err
end


using Optim

od = Optim.OnceDifferentiable(p->errfun(p), p0, autodiff=:finite)
mfit = Optim.optimize(od, lb, ub, p0, Fminbox(LBFGS()))
pmin = Optim.minimizer(mfit)


AbsScale, Tscale = pmin
@show AbsScale
@show Tscale;


using OrdinaryDiffEq
# let grab the derivative of Emax function
import Pumas.IVIVC: e_der
g = e_der
pmin = vitro_data[1]["fast"].pmin
f(c, p, t) = kel * conc_auc * AbsScale * Tscale * g(t * Tscale, pmin) - kel * c
c0 = 0.0
tspan = (vivo_fast_data.time[1], vivo_fast_data.time[end])
prob = ODEProblem(f, c0, tspan)
sol = OrdinaryDiffEq.solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
sol


scatter(vivo_fast_data.time, vivo_fast_data.conc, label="Original C (fast)")
plot!(sol, label="Predicted C (fast)", xlabel="time", ylabel="Plasma Conc. Profile")


cmax_pe, auc_pe = percentage_prediction_error(vivo_fast_data.time, vivo_fast_data.conc,
                              sol.t, sol.u)   # sol.u is just predicted c
@show cmax_pe
@show auc_pe;


model = IVIVCModel(vitro_data, uir_data, vivo_data, vitro_model=:emax, uir_frac=1.0, deconvo_method=:wn, ivivc_model=:three);


ivivc_plot(model)


levy_plot(model)


sol = predict_vivo(model, "fast");


time = vivo_data[1]["fast"].time
conc = vivo_data[1]["fast"].conc
cmax_pe, auc_pe = percentage_prediction_error(time, conc, sol.t, sol.u)
@show cmax_pe, auc_pe;
scatter(time, conc, label="Original C (fast)")
plot!(sol, label="predicted C (fast)", xlabel="time", ylabel="Plasma Conc. Profile")


sol = predict_vivo(model, "medium")
time = vivo_data[1]["medium"].time
conc = vivo_data[1]["medium"].conc
cmax_pe, auc_pe = percentage_prediction_error(time, conc, sol.t, sol.u)
@show cmax_pe, auc_pe;
scatter!(time, conc, label="Original C (medium)")
plot!(sol, label="Predicted C (medium)")


sol = predict_vivo(model, "slow")
time = vivo_data[1]["slow"].time
conc = vivo_data[1]["slow"].conc
cmax_pe, auc_pe = percentage_prediction_error(time, conc, sol.t, sol.u)
@show cmax_pe, auc_pe;
scatter!(time, conc, label="Original C (slow)")
plot!(sol, label="Predicted C (slow)")


to_csv(model, path=homedir())

