
using CSV, DataFrames, PuMaS
using DataFrames: head


nmtran_data = example_nmtran_data("data1")
show(nmtran_data)


data = DataFrame(CSV.File(nmtran_data)) # Make a DataFrame from the text file.
head(data) # Show the first six rows


using PuMaS
@doc process_nmtran


population = process_nmtran(data, [:sex, :wt])
show(population)


using Distributions
distribution = Normal(5.5, 1.5)
show(distribution)


glucose = rand(distribution, 5) # Five draws from the distribution


obs = DataFrame(time = 0:8:32, # Makes a sequence starting at 0 by 8 until 32
                glucose = glucose)
head(obs)


cvs = (sex = "male", age = 25)
show(cvs)


evs = DosageRegimen(5, addl = 3, ii = 8u"hr")
head(evs.data)


subject = Subject(cvs = cvs, obs = obs, evs = evs)
show(subject)


for obs in subject.observations
    println(obs)
end


for evs in subject.events
    println(evs)
end


regimens = DosageRegimen([0.005u"g", 0.0025u"g"], time = [0, 90u"minute"])
head(regimens.data)


subject2 = Subject(id = 2,
                   evs = regimens)
show(subject2)


for evs in subject2.events
    println(evs)
end


@doc Population


ids = 1:6 # There will be 6 subject 1, 2, 3, 4, 5, 6


choose_sex_age() = (sex = rand(["male", "female"]),
                    age = rand(21:25))


cvs = [ choose_sex_age() for i in 1:6 ]
show(cvs)


obs = [ DataFrame(time = 0:4:8, dv = rand(3)) for id in 1:6 ]
head.(obs[1:2])


regimen1 = DosageRegimen(15, addl = 3, ii = 4)
regimen2 = DosageRegimen(10, addl = 4, ii = 3)
regimens = vcat(repeat([regimen1], 3),
                repeat([regimen2], 3))



subjects = Subject[]
for (id, cvs, evs) in zip(ids, cvs, regimens)
    push!(subjects, Subject(id = id, cvs = cvs, evs = evs))
end
show(subjects)


population = Population(subjects)
show(population)


for subject in population
    println(subject)
end

