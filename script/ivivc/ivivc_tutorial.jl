
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

