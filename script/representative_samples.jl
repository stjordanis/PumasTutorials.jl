
using CSV, DataFrames, DataFramesMeta, Distributions, HTTP, InfoZIP, LazyJSON,
      Optim, StatsBase, StatPlots


cross_table =
    string("https://api.census.gov/data/", # Base URL
           "2017/", # Vintage 2017
           "pep/", # Population Estimates
           "charage?", # Characteristics by Single Year of Age
           "get=AGE,HISP,POP,RACE,SEX&", # Variables
           "for=us", # For the whole United States
           ) |>
    (url -> HTTP.request("GET", url)) |> # GET request
    (response -> response.status == 200 ? # Verify status code OK
                 String(response.body) : response.status) # Get body of response


@views function detailedtable(jsontext)
    # We obtain the variable names (suppress the geography)
    vars = Symbol.(collect(LazyJSON.parse(jsontext)[1][1:end - 1]))
    # For each record we obtain a NamedTuple with the values parsed
    parserow(val) = NamedTuple{Tuple(vars)}(parse.(Int, val[1:end - 1]))
    # Function to determine race
    # Latino if Hispanic
    # White if Single Race White (Non-Hispanic)
    # Black if Single Race Black (Non-Hispanic)
    # Asian if Single Race Asian (Non-Hispanic)
    function parserace(hisp, race)
        if hisp == 2
            output = "Latino"
        elseif race == 1
            output = "White"
        elseif race == 2
            output = "Black"
        elseif race == 4
            output = "Asian"
        end
        output
    end
    DataFrame(parserow.(LazyJSON.parse(jsontext)[2:end])) |>
    (df -> @where(df,
                  (:SEX .∈ Ref(1:2)) .& # Male = 1, Female = 2
                  (:AGE .∈ Ref(18:85)) .& # Subset ages 18-85 years old
                  ((:HISP .== 2) .| # Hispanic = 2
                  # Not Hispanic and one of the select races
                  # White = 1, Black = 2, Asian = 4
                   ((:HISP .== 1) .& (:RACE .∈ Ref([1, 2, 4])))))) |>
    (df -> @transform(df,
                      SEX = ifelse.(:SEX .== 1, "Male", "Female") |> categorical,
                      RACE = parserace.(:HISP, :RACE) |> categorical)) |>
    (df -> aggregate(groupby(df, [:SEX, :RACE, :AGE]), sum)) |> # Collapse Latino
    (df -> names!(df[[:SEX, :RACE, :AGE, :POP_sum]], [:sex, :race, :age, :count]))
end


cross_table = detailedtable(cross_table)
head(cross_table)


females = (@linq cross_table |>
           where(:sex .== "Female") |>
           select(:race, :age, :count)) |>
           (df -> unstack(df, :race, :count)) |> # To wide
           (df -> Matrix(df[2:end])) |>
           (A -> A / sum(A)) # conditional joint probabilities
males = (@linq cross_table |>
         where(:sex .== "Male") |> # Filter males
         select(:race, :age, :count)) |>
         (df -> unstack(df, :race, :count)) |> # To wide
         (df -> Matrix(df[2:end])) |>
         (A -> A / sum(A)) # conditional joint probabilities


println("Valid Probabilities: $(sum(females) ≈ sum(males) ≈ 1)")
pepall5n = cat(females, males, dims = 3)
size(pepall5n)


function gen_join_distribution(probabilities, output = Vector{Int}())
    d = ndims(probabilities) # How many dimensions
    iszero(ndims(probabilities)) && return output # If no more dimensions, done
    o = sample(1:size(probabilities, d), # How many options for the next dimensions?
               # Probability weights are given by summing up through the dimensions
               sum.(selectdim.(Ref(probabilities),
                               d,
                               1:size(probabilities, d))) |>
               weights)
    push!(output, o) # Add the observed outcome
    # repeat process given the realization in the previous dimensions
    gen_join_distribution(selectdim(probabilities, d, o), output)
end
labs = [(:sex, levels(cross_table.sex)),
        (:race, levels(cross_table.race)),
        (:age, levels(cross_table.age))]
function gen_sex_race_age(joint_distribution, labs)
    NamedTuple{Tuple(first.(labs))}(getindex.(last.(labs),
                                    gen_join_distribution(joint_distribution)))
end


data = map(x -> gen_sex_race_age(pepall5n, labs), 1:1_000) |>
       DataFrame |>
       categorical!
head(data)


function parserace(race, hisp)
    if hisp ≠ 12
        output = "Latino"
    elseif race == 1
        output = "White"
    elseif race == 2
        output = "Black"
    elseif race == 4
        output = "Asian"
    else
        output = ""
    end
    output
end
nchs = string("ftp://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/",
              "NHIS/",
              "2017/",
              "samadultcsv.zip") |>
    (url -> download(url, joinpath(tempdir(), "NCHS.zip"))) |> # Downloads zip file
    (file -> open_zip(file)["samadult.csv"]) |> # Read the file
    IOBuffer |> # Suitable stream to parse as text file rather than string
    CSV.File |> # Reads coma separated file
    DataFrame |> # Materialize as DataFrame
    (df -> (@linq df |>
            # Exclude currently pregnant and invalid BMI records
            where(ifelse.(ismissing.(:PREGNOW), true, :PREGNOW .== 2) .&
                  (1 .≤ :BMI .≤ 9994)) |>
            select(:SEX, :RACERPI2, :HISPAN_I, :AGE_P, :AHEIGHT, :BMI,
                   :WTFA_SA))) |> # Survey weights
    (df -> @transform(df,
                      SEX = ifelse.(:SEX .== 1, "Male", "Female"),
                      RACERPI2 = parserace.(:RACERPI2, :HISPAN_I),
                      BMI = :BMI / 100)) |> # Adjust implied decimals
    (df -> @select(df, :SEX, :RACERPI2, :AGE_P, :AHEIGHT, :BMI, :WTFA_SA)) |>
    (df -> names!(df, [:sex, :race, :age, :height, :bmi, :wts])) |>
    # Ensure valid values for height given documentation
    # Filter races of interest
    (df -> @where(df,
                  (((:sex .== "Male") .& (59 .≤ :height .≤ 76)) .|
                   ((:sex .== "Female") .& (59 .≤ :height .≤ 70))) .&
                   (:race .≠ ""))) |>
    disallowmissing! |> # Drop missing values
    categorical! # String variables are treated as nominal variables
head(nchs)


function mle_Weibull(data)
    x₀ = fill(2.0, 2) # Initial Parameters
    # Set up the objective function, first and second derivatives
    wts = data.wts ./ sum(data.wts)
    td = TwiceDifferentiable(x -> -sum(logpdf.(Weibull(x[1], # Flip sign
                                                       x[2]),
                                               data.bmi) .*
                                        wts),
                            x₀)
    lx = ones(2) # Lower bound is one for both parameters
    ux = fill(Inf, 2) # No upper bound
    tdc = TwiceDifferentiableConstraints(lx, ux) # Autodiff
    res = optimize(td, tdc, x₀, IPNewton()) # Solve the problem
    Weibull(Optim.minimizer(res)...) # MLE Distribution
end


height_bmi =
    by(nchs, [:sex, :race]) do subdf
        DataFrame(height = fit(Normal, subdf.height, float(subdf.wts)),
                  bmi = mle_Weibull(subdf))
    end
head(height_bmi)


function add_height_weight!(data, height_bmi)
    data.height = 0 # Initialize height
    data.weight = 0 # Initialize weight
    for row in eachrow(data)
        # Unpack
        sex, race, height, bmi =
            row.sex, row.race, row.height, row.bmi
        # Identify relevant cases
        idx = (data.sex .== sex) .& (data.race .== race)
        # Sample height for group
        data.height[idx] = round.(Int, rand(height, count(idx)))
        # Go from BMI to pounds
        data.weight[idx] = round.(Int, rand(bmi, count(idx)) .*
                                       data.height[idx].^2 / 703)
    end
    data
end


add_height_weight!(data, height_bmi) |> head


@views rand(MvNormal([70, 30], [10^2 32.
                                32   5^2]),
            10) |>
    (x -> DataFrame(age = round.(Int, x[1,:]),
                    weight = round.(Int, x[2, :])))

