# PuMaSTutorials.jl

Tutorials for pharmaceutical modeling and simulation with PuMaS.jl

## Current Tutorials

- [PuMaS for Multiple Response Pk/Pd](https://github.com/UMCTM/PuMaSTutorials.jl/blob/master/pdf/multiple_response.pdf)
- [Data Generation for PuMaS](https://github.com/UMCTM/PuMaSTutorials.jl/blob/master/pdf/data_generation.pdf)
- [PBPK in PuMaS, A Model for ACAT](https://github.com/UMCTM/PuMaSTutorials.jl/blob/master/pdf/pbpk_acat.pdf)

# Developer Documentation

## Building the Files

To build the file `multiple_response.jmd`, use the following commands:

```julia
using PuMaSTutorials
PuMaSTutorials.weave_file("multiple_response.jmd")
```

To build all of the files, do `PuMaSTutorials.weave_all()`.
