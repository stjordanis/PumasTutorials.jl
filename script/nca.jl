
using PuMaS.NCA


using PuMaS, CSV

file = PuMaS.example_nmtran_data("nca_test_data/dapa_IV")
data = CSV.read(file)


first(data, 6) # take first 6 rows


NCA.auc(data[:CObs][1:16], data[:TIME][1:16])


NCA.auc(data[:CObs][1:16], data[:TIME][1:16], method=:linuplogdown)


pop = parse_ncadata(data, id=:ID, time=:TIME, conc=:CObs, amt=:AMT_IV, formulation=:Formulation, iv="IV", llq=0)
pop[1]


NCA.auc(pop)


NCA.auc(pop[2], auctype=:last)


NCA.auc(pop, auctype=:last)


NCA.auc(pop[1], interval=(10,Inf))


NCA.auc(pop[1], interval=[(10,Inf), (10, 15)])


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


NCA.lambdaz(pop[1], slopetimes=[1,2,3])


NCA.tmax(pop[1])
NCA.cmax(pop[1])
NCA.cmax(pop[1], interval=(20, 24))
NCA.cmax(pop[1], interval=[(20, 24), (10, 15)])


NCA.tlast(pop[1])
NCA.clast(pop[1])


NCA.thalf(pop[1])


NCA.interpextrapconc(pop[1], 12., interpmethod=:linear)

