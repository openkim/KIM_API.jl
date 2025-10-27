# Getting Started

This guide will help you get up and running with KIM_API.jl.

## Installation

### Step 1: Install KIM-API

KIM-API provides packages for various OS and distributions. For all the information, and instructions on how to install KIM-API, visit [KIM-API Github repo](https://github.com/openkim/kim-api), and the [installation notes](https://github.com/openkim/kim-api/blob/master/INSTALL).


Easiest way to install KIM-API is via Conda:

```bash
conda create -n kim-api kim-api=2.4 -c conda-forge
conda activate kim-api
export KIM_API_LIB=${CONDA_PREFIX}/lib/libkim-api.so # for KIM_API.jl to find the library
```

### Step 2: Install KIM_API.jl

```julia
using Pkg
Pkg.add(url="https://github.com/openkim/KIM_API.jl.git")
```

### Step 3: Test Installation

```julia
using KIM_API

# This should not throw an error
println("KIM-API library loaded successfully!")
```

If you get a library loading error, you may need to set the `KIM_API_LIB` environment variable:

```bash
export KIM_API_LIB="/path/to/libkim-api.so"
```

## Your First Calculation

Let's compute the energy and forces for a simple silicon system using the Stillinger-Weber potential.

```julia
using KIM_API
using StaticArrays
using LinearAlgebra

# Create a KIM model
model = KIM_API.KIMModel("SW_StillingerWeber_1985_Si__MO_405512056662_006")

# Define a simple silicon system
species = ["Si", "Si"]
positions = [
        SVector(0.    , 0.    , 0.    ),
        SVector(1.3575, 1.3575, 1.3575),
]
# Silicon lattice cell
cell = Matrix([[0.0 2.715 2.715] 
               [2.715 0.0 2.715] 
               [2.715 2.715 0.0]])

# Use periodic boundary conditions
pbc = [true, true, true]

# Compute energy and forces
results = model(species, positions, cell, pbc)

println("Energy: $(results[:energy]) eV")
println("Forces:")
for (i, force) in enumerate(eachcol(results[:forces]))
    println("  Atom $i: $force eV/Å")
end
```

Expected output:
```
Energy: -8.67279650983989 eV

Forces:
  Atom 1: [0.0, 0.0, -4.336808689942018e-19] eV/Å
  Atom 2: [-6.505213034913027e-19, -6.505213034913027e-19, -8.673617379884035e-19] eV/Å
```

## Understanding the Results

- **Energy**: Total potential energy of the system in eV
- **Forces**: Forces on each atom in eV/Å (note the units depend on your chosen unit system)

The forces should sum to approximately zero for an isolated system due to Newton's third law.

## Exploring Available Models

KIM-API provides access to hundreds of validated interatomic models. You can browse them at [openkim.org](https://openkim.org).

### Working with Different Elements

```julia
# Copper using EAM potential
cu_model = KIM_API.KIMModel("EAM_Dynamo_MendelevSordeletKramer_2007_CuZr__MO_120596890176_005")

# Define copper system
cu_species = ["Cu", "Cu", "Cu", "Cu"]
cu_positions = [
    SVector(0.0, 0.0, 0.0),
    SVector(1.805, 1.805, 0.0),
    SVector(1.805, 0.0, 1.805),
    SVector(0.0, 1.805, 1.805)
]

# FCC copper lattice parameter
a_cu = 3.61
cu_cell = Matrix(a_cu * I(3))
cu_pbc = [true, true, true]

cu_results = cu_model(cu_species, cu_positions, cu_cell, cu_pbc)
println("Copper energy: $(cu_results[:energy]) eV")
```
```
Copper energy: -13.120417106582778 eV
```

## Unit Systems

KIM_API.jl supports multiple unit systems. The default is `:metal` (LAMMPS metal units):

```julia
# Explicit unit specification
model_real = KIM_API.KIMModel("SW_StillingerWeber_1985_Si__MO_405512056662_006",
                      units=:real)  # kcal/mol, Å, fs

model_si = KIM_API.KIMModel("SW_StillingerWeber_1985_Si__MO_405512056662_006",
                    units=:si)    # J, m, s
```

| Unit System | Length | Energy | Time | Notes |
|-------------|--------|--------|------|-------|
| `:metal`    | Å      | eV     | ps   | Most common for atomistic simulations |
| `:real`     | Å      | kcal/mol | fs | LAMMPS real units |
| `:si`       | m      | J      | s    | SI base units |
| `:cgs`      | cm     | erg    | s    | CGS units |
| `:electron` | Bohr   | Hartree| fs   | Atomic units |

You can also specify custom collection of units by passing a tuple of units during model creation:

```julia
custom_units = (length=KIM_API.A, energy=KIM_API.eV, time=KIM_API.fs, charge=KIM_API.e, temperature=KIM_API.K)
model_custom = KIM_API.KIMModel("SW_StillingerWeber_1985_Si__MO_405512056662_006",
                        units=custom_units)
```

## Performance Considerations

### Reuse Model Functions

Creating a KIM model has overhead, so reuse the function:

```julia
# Good: create once
model = KIM_API.KIMModel("SW_StillingerWeber_1985_Si__MO_405512056662_006")

# Use many times
for configuration in configurations
    results = model(species, configuration.positions, cell, pbc)
    # Process results...
end
```

### Choose Appropriate Compute Options

If you only need energy (not forces):

```julia
model_energy_only = KIM_API.KIMModel("SW_StillingerWeber_1985_Si__MO_405512056662_006",
                             compute=[:energy])
```

This can be significantly faster for large systems.

## Common Issues

### Library Not Found

```
ERROR: KIM-API library not found in system paths.
```

**Solution**: Install KIM-API or set the environment variable:
```bash
export KIM_API_LIB="/usr/local/lib/libkim-api.so"
```

### Model Not Found

```
ERROR: Model creation failed with error code: 1
```

**Solutions**:
1. Check the model name spelling
2. Install the model: `kim-api-collections-management install user <model-name>`
3. Browse available models at [openkim.org](https://openkim.org)

### Species Not Supported

```
ERROR: Species 'Unobtainium' not supported by model
```

**Solution**: Check which elements the model supports. Most models specify supported elements in their documentation.

## Next Steps

- Check out advanced usage patterns in the examples
- Check out [Examples](examples.md) for real-world applications
- Explore the [API Reference](api/highlevel.md) for detailed documentation
- Visit [openkim.org](https://openkim.org) to discover more models
