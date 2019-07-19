
using Pumas.NCA


using Pumas, CSV

file = Pumas.example_nmtran_data("nca_test_data/dapa_IV")
data = CSV.read(file);


first(data, 6) # take first 6 rows


NCA.auc(data[:CObs][1:16], data[:TIME][1:16])


NCA.auc(data[:CObs][1:16], data[:TIME][1:16], method=:linuplogdown)


timeu = u"hr"
concu = u"mg/L"
amtu  = u"mg"
data.id = data.ID
data.time = data.TIME
data.conc = data.CObs
data.amt = data.AMT_IV
data.route = "iv"
pop = read_nca(data, llq=0concu, timeu=timeu, concu=concu, amtu=amtu)


NCA.auc(pop)


NCA.auc(pop[2], auctype=:last)


NCA.auc(pop, auctype=:last)


NCA.auc(pop[1], interval=(10,Inf).*timeu)


NCA.auc(pop[1], interval=[(10,Inf).*timeu, (10, 15).*timeu])


NCA.auc_extrap_percent(pop[1])


NCA.aumc_extrap_percent(pop[1])
NCA.aumc(pop[1])


NCA.lambdaz(pop[1])


NCA.lambdazr2(pop)
NCA.lambdazadjr2(pop)
NCA.lambdazintercept(pop)
NCA.lambdaztimefirst(pop)
NCA.lambdaznpoints(pop)


NCA.lambdaz(pop[1], threshold=3)


NCA.lambdaz(pop[1], idxs=[10, 15, 16])


NCA.lambdaz(pop[1], slopetimes=[1,2,3].*timeu)


NCA.tmax(pop[1])
NCA.cmax(pop[1])
NCA.cmax(pop[1], interval=(20, 24).*timeu)
NCA.cmax(pop[1], interval=[(20, 24).*timeu, (10, 15).*timeu])


NCA.tlast(pop[1])
NCA.clast(pop[1])


NCA.thalf(pop[1])


NCA.interpextrapconc(pop[1], 12timeu, method=:linear)


using Plots # load the plotting library
plot(pop)


plot(pop, loglinear=false)


plot(pop, linear=false)


report = NCAReport(pop)


reportdf = NCA.to_dataframe(report)
first(reportdf,6) # Print only the first 6 rows


multiple_doses_file = Pumas.example_nmtran_data("nca_test_data/dapa_IV_ORAL")
mdata = CSV.read(multiple_doses_file)

timeu = u"hr"
concu = u"mg/L"
amtu  = u"mg"
mdata.id = mdata.ID
mdata.time = mdata.TIME
mdata.conc = mdata.COBS
mdata.amt = mdata.AMT
mdata.route = replace(mdata.FORMULATION, "IV"=>"iv", "ORAL"=>"ev")
mdata.occasion = mdata.OCC
mpop = read_nca(mdata, timeu=timeu, concu=concu, amtu=amtu)


plot(mpop)


NCA.auc(mpop)


rep = NCAReport(mpop, ithdose=1)
reportdf = NCA.to_dataframe(rep)
first(reportdf,6) # Print only the first 6 rows

