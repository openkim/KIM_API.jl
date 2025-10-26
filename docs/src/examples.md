# Examples

This page provides practical examples of using KIM_API.jl for various molecular simulation tasks.

## Basic Examples

### Computing Energy and Forces with KIMModel

```julia
using KIM_API
using StaticArrays
using LinearAlgebra

# Create a KIM model for silicon
model = KIMModel("SW_StillingerWeber_1985_Si__MO_405512056662_006")

# Define a silicon dimer
species = ["Si", "Si"]
positions = [
    SVector(0.0, 0.0, 0.0),     # First silicon atom
    SVector(2.35, 0.0, 0.0)     # Second silicon atom at equilibrium distance
]

# Set up periodic cell (large enough to avoid self-interaction)
cell = Matrix(10.0 * I(3))  # 10×10×10 Å cubic cell
pbc = [true, true, true]

# Compute properties
results = model(species, positions, cell, pbc)

println("Energy: $(results[:energy]) eV")
println("Force on atom 1: $(results[:forces][:, 1]) eV/Å")
println("Force on atom 2: $(results[:forces][:, 2]) eV/Å")
```

### Using AtomsBase and AtomsCalculators

```julia
using KIM_API, AtomsBase, AtomsCalculators
using StaticArrays, Unitful

# Create a KIM calculator
calc = KIMCalculator("SW_StillingerWeber_1985_Si__MO_405512056662_006")

# Create AtomsBase system with units
particles = [
    :Si => SVector(0.0u"Å", 0.0u"Å", 0.0u"Å"),
    :Si => SVector(2.35u"Å", 0.0u"Å", 0.0u"Å")
]

cell_vectors = (
    SVector(10.0u"Å", 0.0u"Å", 0.0u"Å"),
    SVector(0.0u"Å", 10.0u"Å", 0.0u"Å"),
    SVector(0.0u"Å", 0.0u"Å", 10.0u"Å")
)

system = FlexibleSystem(particles; cell_vectors=cell_vectors, periodicity=(true, true, true));

# Use AtomsCalculators interface
energy = AtomsCalculators.potential_energy(calc, system)
forces = AtomsCalculators.forces(calc, system)

println("Energy: $(energy) eV")
println("Forces: $(forces)")

# Or call calculator directly for all properties
results = calc(system)
println("All results: $(results)")
```

### Energy Minimization

```julia
using KIM_API, Optim
using StaticArrays, LinearAlgebra

# Set up model and initial structure
model = KIM_API.KIMModel("SW_StillingerWeber_1985_Si__MO_405512056662_006")
species = ["Si", "Si"]
cell = [5.43 0.0 0.0; 0.0 5.43 0.0; 0.0 0.0 5.43]
pbc = [true, true, true]

# Objective function for optimization
function objective(x)
    pos = [SVector(0.0, 0.0, 0.0), SVector(x[1], x[2], x[3])]
    results = model(species, pos, cell, pbc)
    return results[:energy]
end

# Optimize second atom position
initial_pos = [2.5, 0.0, 0.0]  # Starting guess
result = result = optimize(objective, initial_pos, BFGS(), Optim.Options(f_tol=1e-4))

println("Optimized position: ", result.minimizer)
println("Minimum energy: ", result.minimum)
```

## Advanced Examples

### Equivalence Between Methods

This example demonstrates that `KIMModel` (raw arrays) and `KIMCalculator` (AtomsBase) produce identical results:

```julia
using KIM_API, AtomsBase, StaticArrays

model_name = "SW_StillingerWeber_1985_Si__MO_405512056662_006"

# Method 1: Raw arrays with KIMModel
model = KIMModel(model_name)
species = ["Si", "Si"]
positions = [SVector(0.0, 0.0, 0.0), SVector(2.35, 0.0, 0.0)]
cell = [5.43 0.0 0.0; 0.0 5.43 0.0; 0.0 0.0 5.43]
pbc = [true, true, true]

raw_results = model(species, positions, cell, pbc)

# Method 2: AtomsBase system with KIMCalculator
calc = KIMCalculator(model_name)
particles = [
    :Si => SVector(0.0u"Å", 0.0u"Å", 0.0u"Å"),
    :Si => SVector(2.35u"Å", 0.0u"Å", 0.0u"Å")
]
cell_vectors = (
    SVector(5.43u"Å", 0.0u"Å", 0.0u"Å"),
    SVector(0.0u"Å", 5.43u"Å", 0.0u"Å"),
    SVector(0.0u"Å", 0.0u"Å", 5.43u"Å")
)
system = FlexibleSystem(particles; cell_vectors=cell_vectors, periodicity=(true, true, true));
atomsbase_results = calc(system)

# Verify equivalence
@assert raw_results[:energy] ≈ atomsbase_results[:energy]
@assert raw_results[:forces] ≈ atomsbase_results[:forces]

println("Both methods produce identical results!")
println("Energy: $(raw_results[:energy]) eV")
```

### Different Unit Systems

```julia
using KIM_API, StaticArrays, LinearAlgebra

model_name = "SW_StillingerWeber_1985_Si__MO_405512056662_006"
species = ["Si", "Si"]
positions = [SVector(0.0, 0.0, 0.0), SVector(2.35, 0.0, 0.0)]
cell = Matrix(5.43 * I(3))
pbc = [true, true, true]

# Different unit systems
units_list = [:metal, :real, :si]
results = Dict()

for units in units_list
    calc = KIMCalculator(model_name, units=units)
    model_fn = calc.model_fn
    result = model_fn(species, positions, cell, pbc)
    results[units] = result[:energy]
end

println("Energy in different units:")
println("  Metal (eV): $(results[:metal])")
println("  Real (kcal/mol): $(results[:real])")
println("  SI (J): $(results[:si])")

# Convert real units to eV for comparison
energy_real_in_eV = results[:real] * 0.043364
println("  Real-> eV: $(energy_real_in_eV)")
```


