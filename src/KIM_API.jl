"""
    KIM_API.jl

A comprehensive Julia interface to the KIM-API (Knowledgebase of Interatomic Models).

This package provides both low-level and high-level interfaces to KIM-API,
enabling Julia users to access the extensive collection of validated 
interatomic models available through the OpenKIM framework.

# Features
- High-level functional interface with `KIMModel()`
- Automatic species mapping and validation
- Efficient neighbor list generation with periodic boundary conditions
- Support for multiple unit systems (metal, real, SI, CGS, atomic)
- Memory-safe wrappers around KIM-API C functions
- Comprehensive error handling and validation

# Quick Start
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

# Modules
- `libkim.jl`: KIM-API library loading and initialization
- `constants.jl`: KIM-API constants, enumerations, and unit systems
- `model.jl`: Low-level model creation and management
- `species.jl`: Species handling and validation
- `neighborlist.jl`: Neighbor list generation and callbacks
- `highlevel.jl`: High-level user interface

# Units
The package supports multiple unit systems:
- `:metal`: Ångström, eV, electron charge, Kelvin, picosecond (LAMMPS metal)
- `:real`: Ångström, kcal/mol, electron charge, Kelvin, femtosecond (LAMMPS real)
- `:si`: meter, Joule, Coulomb, Kelvin, second
- `:cgs`: centimeter, erg, statCoulomb, Kelvin, second
- `:electron`: Bohr, Hartree, electron charge, Kelvin, femtosecond

# KIM-API Integration
This package wraps the KIM-API C library and provides:
- Automatic memory management
- Julia-native data structures
- Type-safe function interfaces
- Comprehensive documentation

# Dependencies
- KIM-API C library (must be installed separately)
- KIMNeighborList.jl for efficient C++ neighbor searching
- StaticArrays.jl for position vectors
- CEnum.jl for C enumeration bindings

For more information about KIM-API and available models, visit:
- https://openkim.org
- https://kim-api.readthedocs.io
"""
module KIM_API
using CEnum
using StaticArrays
using KIMNeighborList
using LinearAlgebra
using Libdl

include("utils.jl") # Utility functions for scatter-add and reduce operations

include("libkim.jl") # Load the KIM API library

include("constants.jl") # Default enumerations in kim-api

include("model.jl") # functions for model creation and destruction
# also includes ComputeArguments stuff for now

include("species.jl") # All the species related functions

include("neighborlist.jl") # Functions for creating neighbor lists

include("highlevel.jl") # Actual high-level API for model computation

function __init__()
    refresh_libkim!(warn_if_missing = true)
end
end
