using Test
using StaticArrays
using LinearAlgebra
using Unitful
using AtomsBase: FlexibleSystem
using AtomsCalculators

@testset "KIMCalculator Integration Tests" begin
    model_name = "SW_StillingerWeber_1985_Si__MO_405512056662_006"

    try
        calc = KIM_API.KIMCalculator(model_name)

        species = ["Si", "Si"]
        positions = [SVector(0.0, 0.0, 0.0), SVector(2.35, 0.0, 0.0)]
        cell_length = 10.0
        cell = Matrix(cell_length * I(3))
        pbc = [true, true, true]

        results = calc(species, positions, cell, pbc)
        @test results isa NamedTuple
        @test isfinite(results[:energy])
        @test size(results[:forces]) == (3, length(species))
        @test all(isfinite, results[:forces])

        @test_throws ErrorException calc(species, positions, cell, pbc; unsupported = 1)

        particles = [
            :Si => SVector(0.0u"Å", 0.0u"Å", 0.0u"Å"),
            :Si => SVector(2.35u"Å", 0.0u"Å", 0.0u"Å"),
        ]
        cell_vectors = (
            SVector(cell_length * u"Å", 0.0u"Å", 0.0u"Å"),
            SVector(0.0u"Å", cell_length * u"Å", 0.0u"Å"),
            SVector(0.0u"Å", 0.0u"Å", cell_length * u"Å"),
        )
        system = FlexibleSystem(particles; cell_vectors = cell_vectors, periodicity = (true, true, true))

        system_results = calc(system)
        @test system_results isa NamedTuple
        @test isapprox(system_results[:energy], results[:energy]; atol = 1e-8)
        @test isapprox(system_results[:forces], results[:forces]; atol = 1e-8)

        energy = AtomsCalculators.potential_energy(calc, system)
        forces = AtomsCalculators.forces(calc, system)
        @test energy ≈ system_results[:energy]
        @test isapprox(forces, system_results[:forces]; atol = 1e-8)

        @testset "Molly Integration" begin
            molly_available = try
                @eval import Molly
                true
            catch
                false
            end

            if molly_available
                ext = Base.get_extension(KIM_API, :KIM_APIMollyExt)
                @test ext !== nothing

                has_forces! = any(
                    m -> m.sig <:
                        Tuple{typeof(AtomsCalculators.forces!), Any, Molly.System, KIM_API.KIMCalculator},
                    methods(AtomsCalculators.forces!),
                )
                @test has_forces!

                has_energy = any(
                    m -> m.sig <:
                        Tuple{typeof(AtomsCalculators.potential_energy), Molly.System, KIM_API.KIMCalculator},
                    methods(AtomsCalculators.potential_energy),
                )
                @test has_energy
            else
                @test_skip "Molly not available in current environment"
            end
        end
    catch e
        if occursin("Model creation failed", sprint(showerror, e)) ||
           occursin("not found", sprint(showerror, e))
            @test_skip "Test model $model_name not available: $e"
        else
            rethrow(e)
        end
    end
end
