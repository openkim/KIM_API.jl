module KIM_APIMollyExt

using KIM_API
import AtomsCalculators
import Molly
using StaticArrays: SVector
using LinearAlgebra: Diagonal
using Unitful: Å, ustrip

@inline function _kim_inputs(sys::Molly.System)
    natoms = length(sys.atoms)
    species = Vector{String}(undef, natoms)
    positions = Vector{SVector{3, Float64}}(undef, natoms)

    @inbounds for i in 1:natoms
        species[i] = String(getfield(sys.atoms[i], :atom_type))
        positions[i] = SVector{3, Float64}(ustrip.(Å, sys.coords[i]))
    end

    sides = Float64.(ustrip.(Å, sys.boundary.side_lengths))
    cell = Matrix(Diagonal(sides))
    pbc = Vector{Bool}(undef, length(sides))
    @inbounds for i in eachindex(sides)
        pbc[i] = !isinf(sides[i])
    end

    return (; species, positions, cell, pbc)
end

function AtomsCalculators.forces(
    sys::Molly.System,
    calc::KIM_API.KIMCalculator;
    kwargs...,
)
    inputs = _kim_inputs(sys)
    results = calc(inputs.species, inputs.positions, inputs.cell, inputs.pbc)
    fu = sys.force_units
    forces = results[:forces]
    n_atoms = size(forces, 2)
    out = Vector{typeof(SVector{3, Float64}(0, 0, 0) * fu)}(undef, n_atoms)

    @inbounds for (i, col) in enumerate(eachcol(forces))
        out[i] = SVector{3, Float64}(col) * fu
    end

    return out
end

function AtomsCalculators.forces!(
    fs,
    sys::Molly.System,
    calc::KIM_API.KIMCalculator;
    kwargs...,
)
    inputs = _kim_inputs(sys)
    results = calc(inputs.species, inputs.positions, inputs.cell, inputs.pbc)
    fu = sys.force_units
    forces = results[:forces]

    @inbounds for (i, col) in enumerate(eachcol(forces))
        fs[i] += SVector{3, Float64}(col) * fu
    end
    return fs
end

function AtomsCalculators.potential_energy(
    sys::Molly.System,
    calc::KIM_API.KIMCalculator;
    kwargs...,
)
    inputs = _kim_inputs(sys)
    energy = calc(inputs.species, inputs.positions, inputs.cell, inputs.pbc)[:energy]
    return energy * sys.energy_units
end

end # module
