# PumasTutorials.jl

PumasTutorials.jl holds PDFs, webpages, and interactive Jupyter notebooks
showing how to do pharmaceutical modeling and simulation with Pumas.jl.

## Interactive Notebooks

To run the tutorials interactively via Jupyter notebooks, install the package
and open the tutorials like:

```julia
using Pkg
pkg"add https://github.com/UMCTM/PumasTutorials.jl"
using PumasTutorials
PumasTutorials.open_notebooks()
```

## Table of Contents

- Introduction
  - [Introduction to Pumas](https://github.com/UMCTM/PumasTutorials.jl/blob/master/pdf/introduction/introduction.pdf)
  - [Generating and Simulating Populations](https://github.com/UMCTM/PumasTutorials.jl/blob/master/pdf/introduction/simulating_populations.pdf)
  
  - [Introduction to Noncompartmental Analysis (NCA)](https://github.com/UMCTM/PumasTutorials.jl/blob/master/pdf/nca/basic_nca.pdf)
- Models
  - [PBPK in Pumas, A Model for ACAT](https://github.com/UMCTM/PumasTutorials.jl/blob/master/pdf/pbpk/pbpk_acat.pdf)

# Developer Documentation

## Contributing

First of all, make sure that your current directory is `PumasTutorials`. All
of the files are generated from the Weave.jl files in the `tutorials` folder.
To run the generation process, do for example:

```julia
using Pkg, PumasTutorials
cd(joinpath(dirname(pathof(PumasTutorials)), ".."))
Pkg.pkg"activate ."
Pkg.pkg"instantiate"
PumasTutorials.weave_file("introduction","introduction.jmd")
```

To generate all of the notebooks, do:

```julia
PumasTutorials.weave_all()
```

If you add new tutorials which require new packages, simply updating your local
environment will change the project and manifest files. When this occurs, the
updated environment files should be included in the PR.
