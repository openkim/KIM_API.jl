# KIM_API.jl

<p align="center">
<img src="./kimapijl.png" alt="KIM API JL Logo" width="300" />
</p>

Julia interface to the [KIM-API](https:https://kim-api.readthedocs.io) (Knowledgebase of Interatomic Models). 
This is a low-level interface to the KIM-API, allowing you to access interatomic models directly from Julia.
Think of it as the Julia equivalent of the KIMPY Python package.

[Documentation](https://openkim.github.io/KIM_API.jl/)

## Installation

For latest version:
```julia
using Pkg
Pkg.add(url="https://github.com/openkim/KIM_API.jl")
```

For stable version:
```julia
using Pkg
Pkg.add("KIM_API")
```

> This package was earlier called KIMPortableModels.jl

## Quick Start

Export the location of the KIM-API library:

```shell
export KIM_API_LIB=/path/to/libkim-api.so
```

Then, you can use the package as follows:

```julia
using KIM_API, StaticArrays, LinearAlgebra

# Create model function
model = KIM_API.KIMModel("SW_StillingerWeber_1985_Si__MO_405512056662_006")

# Define system
species = ["Si", "Si"]
positions = [
    SVector(0.    , 0.    , 0.    ),
    SVector(1.3575, 1.3575, 1.3575),
]
cell = Matrix([[0.0 2.715 2.715] 
               [2.715 0.0 2.715] 
               [2.715 2.715 0.0]])
pbc = [true, true, true]

# Compute properties
results = model(species, positions, cell, pbc)
println("Energy: ", results[:energy])
println("Forces: ", results[:forces])
```

## Integration with Molly.jl

You can directly use `KIM_API` calculators as general interactions in Molly.jl simulations:

```julia
using Molly, KIM_API, StaticArrays, Unitful, UnitfulAtomic

calc = KIM_API.KIMCalculator("SW_StillingerWeber_1985_Si__MO_405512056662_006";
                              units=:metal)
sys = System(atoms = fill(Atom(atom_type="Si", mass=28.0855u"u"), 2),
             coords = [SVector(0.,0.,0.), SVector(3.,0.,0.)] .* u"Å",
             boundary = CubicBoundary(20.0u"Å"),
             general_inters = (kim = calc,),
             force_units = u"eV/Å", 
             energy_units = u"eV")

println(forces(sys), potential_energy(sys))
```

## Features

- Access to all KIM models
- Automatic neighbor list generation
- Support for periodic boundary conditions
- Multiple unit systems (metal, real, SI, CGS, electron)

## Requirements

- Julia 1.10+
- KIM-API library (for model calculations)
- KIMNeighborList.jl (C++ backend), StaticArrays.jl


## Documentation

Full documentation is available at [https://openkim.github.io/KIM_API.jl/](https://openkim.github.io/KIM_API.jl/)

## Testing

Run the test suite with:

```julia
using Pkg
Pkg.test("KIM_API")
```

## TODO
- Move to 1 based numbering internally for consistency 
- Performance optimizations
- Additional model features support
- Test ML models

## License

MIT
