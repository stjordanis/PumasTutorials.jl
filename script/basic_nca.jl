
using PuMaS.NCA


using PuMaS, CSV

file = PuMaS.example_nmtran_data("nca_test_data/dapa_IV")
data = CSV.read(file)


first(data, 6) # take first 6 rows


NCA.auc(data[:CObs][1:16], data[:TIME][1:16])


NCA.auc(data[:CObs][1:16], data[:TIME][1:16], method=:linuplogdown)


timeu = u"hr"
concu = u"mg/L"
amtu  = u"mg"
pop = parse_ncadata(data, id=:ID, time=:TIME, conc=:CObs, amt=:AMT_IV, formulation=:Formulation, iv="IV",
  llq=0concu, timeu=timeu, concu=concu, amtu=amtu)


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


NCA.lambdaz(pop[1], threshold=2)


NCA.lambdaz(pop[1], idxs=[10, 15, 16])


NCA.lambdaz(pop[1], slopetimes=[1,2,3].*timeu)


NCA.tmax(pop[1])
NCA.cmax(pop[1])
NCA.cmax(pop[1], interval=(20, 24).*timeu)
NCA.cmax(pop[1], interval=[(20, 24).*timeu, (10, 15).*timeu])


NCA.tlast(pop[1])
NCA.clast(pop[1])


NCA.thalf(pop[1])


NCA.interpextrapconc(pop[1], 12timeu, interpmethod=:linear)


using Plots # load the plotting library
plot(pop)


plot(pop, loglinear=false)


plot(pop, linear=false)


report = NCAReport(pop)


NCA.to_dataframe(report)


multiple_doses_file = PuMaS.example_nmtran_data("nca_test_data/dapa_IV_ORAL")
mdata = CSV.read(multiple_doses_file)

timeu = u"hr"
concu = u"mg/L"
amtu  = u"mg"
mpop = parse_ncadata(mdata, time=:TIME, conc=:COBS, amt=:AMT, formulation=:FORMULATION, occasion=:OCC,
                                     iv="IV", timeu=timeu, concu=concu, amtu=amtu)


plot(mpop)


NCA.auc(mpop)


rep = NCAReport(mpop, ithdose=1)
NCA.to_dataframe(rep)

