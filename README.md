# PuMaSTutorials.jl

Tutorials for pharmaceutical modeling and simulation with PuMaS.jl

## Current Tutorials

### Introduction

- [Introduction to PuMaS](https://github.com/UMCTM/PuMaSTutorials.jl/blob/master/pdf/introduction/multiple_response.pdf)
- [Generating and Simulating Populations](https://github.com/UMCTM/PuMaSTutorials.jl/blob/master/pdf/introduction/data_generation.pdf)
- [Introduction to Noncompartmental Analysis (NCA)](https://github.com/UMCTM/PuMaSTutorials.jl/blob/master/pdf/nca/basic_nca.pdf)

### Models

- [PBPK in PuMaS, A Model for ACAT](https://github.com/UMCTM/PuMaSTutorials.jl/blob/master/pdf/pbpk/pbpk_acat.pdf)

# Developer Documentation

## Building the Files

To build the file `multiple_response.jmd`, use the following commands:

```julia
using PuMaSTutorials
PuMaSTutorials.weave_file("multiple_response.jmd")
```

To build all of the files, do `PuMaSTutorials.weave_all()`.
