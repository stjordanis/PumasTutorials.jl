
using Bioequivalence


using CSV, DataFrames, StatsBase, Bioequivalence


PJ31 = CSV.read(string(dirname(pathof(Bioequivalence)),
                       "/../data/Nonparametric/PJ2006_3_1.tsv"))
first(PJ31, 6)


Tmax = read_be(PJ31, :Tmax)


Tmax.model


PJ46 = CSV.read(string(dirname(pathof(Bioequivalence)),
                       "/../data/Williams/PJ2006_4_6.tsv"))
first(PJ46, 6)


NP = read_be(PJ46, nonparametric = true)


ClaytonandLeslie1981 = CSV.read(string(dirname(pathof(Bioequivalence)),
                                       "/../data/Parallel/FSL2015_1.tsv"))
first(ClaytonandLeslie1981, 6)


# notice it defaults to the AUC endpoint
Parallel = read_be(ClaytonandLeslie1981) # Same as read_be(ClaytonandLeslie1981, :AUC)


Parallel.model


Parallel.result


coeftable(Parallel)


SLF2014 = CSV.read(string(dirname(pathof(Bioequivalence)),
                          "/../data/2S2P/SLF2014_8.tsv"))
first(SLF2014, 6)


# notice it defaults to the AUC endpoint
Crossover = read_be(SLF2014)


loglikelihood(Crossover.model.model)


Crossover.result


ChowandLiu2009 = CSV.read(string(dirname(pathof(Bioequivalence)),
                                 "/../data/Balaam/CL2009_9_2_1.tsv"))
first(ChowandLiu2009, 6)


# notice it defaults to the AUC endpoint
Balaam = read_be(ChowandLiu2009)


Balaam.model.σ


PJ41 = CSV.read(string(dirname(pathof(Bioequivalence)),
                       "/../data/Dual/PJ2006_4_1.tsv"))
first(PJ41, 6)


Dual = read_be(PJ41, :Cmax)


PJ43 = CSV.read(string(dirname(pathof(Bioequivalence)),
                       "/../data/2S4P/PJ2006_4_3.tsv"))
first(PJ43, 6)


Inner = read_be(PJ43)


PJ44 = CSV.read(string(dirname(pathof(Bioequivalence)),
                       "/../data/2S4P/PJ2006_4_4.tsv"))
first(PJ44, 6)


Outer = read_be(PJ44, :Cmax)


PJ45 = CSV.read(string(dirname(pathof(Bioequivalence)),
                       "/../data/Williams/PJ2006_4_5.tsv"))
first(PJ45, 6)


W3F = read_be(PJ45)


W3F = read_be(PJ45, reference = 'S')


PJ46 = CSV.read(string(dirname(pathof(Bioequivalence)),
                       "/../data/Williams/PJ2006_4_6.tsv"))
first(PJ46, 6)


W4F = read_be(PJ46)


@doc BioequivalenceStudy


@doc generate_design

