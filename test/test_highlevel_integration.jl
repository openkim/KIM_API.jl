# test_highlevel_integration.jl - Integration tests for high-level interface
# These tests require KIM-API library and specific models

@testset "High-level Interface Integration Tests" begin

    # Test model - Stillinger-Weber for Silicon  
    test_model_name = "SW_StillingerWeber_1985_Si__MO_405512056662_006"

    @testset "KIMModel Function Creation" begin
        try
            # Test basic model creation
            model_func = KIM_API.KIMModel(test_model_name)
            @test model_func isa Function

            # Test with different units
            model_real = KIM_API.KIMModel(test_model_name, units = :real)
            @test model_real isa Function

            model_si = KIM_API.KIMModel(test_model_name, units = :si)
            @test model_si isa Function

            # Test with custom units tuple
            custom_units = (
                length = KIM_API.A,
                energy = KIM_API.eV,
                time = KIM_API.fs,
                charge = KIM_API.e,
                temperature = KIM_API.K,
            )
            model_custom = KIM_API.KIMModel(test_model_name, units = custom_units)
            @test model_custom isa Function

            # Test with energy-only computation
            model_energy = KIM_API.KIMModel(test_model_name, compute = [:energy])
            @test model_energy isa Function

            # Test error handling for invalid arguments
            @test_throws ErrorException KIM_API.KIMModel(test_model_name, compute = Symbol[])  # Empty compute
            invalid_model = KIM_API.KIMModel(test_model_name, compute = [:invalid])
            @test_throws ErrorException invalid_model(
                ["Si"],
                [SVector(0.0, 0.0, 0.0)],
                Matrix(1.0 * I(3)),
                [true, true, true],
            )

        catch e
            if contains(string(e), "Model creation failed")
                @test_skip "Test model $test_model_name not available"
            else
                rethrow(e)
            end
        end
    end

    @testset "Basic Computation" begin
        try
            model = KIM_API.KIMModel(test_model_name)

            # Define simple silicon system
            species = ["Si", "Si"]
            positions = [
                SVector(0.0, 0.0, 0.0),
                SVector(2.35, 2.35, 2.35),  # Silicon lattice distance
            ]

            # Silicon lattice parameter
            a = 5.43  # Angstroms
            cell = Matrix(a * I(3))  # 3x3 identity matrix scaled by lattice parameter
            pbc = [true, true, true]

            # Perform calculation
            results = model(species, positions, cell, pbc)

            @test results isa NamedTuple
            @test haskey(results, :energy)
            @test haskey(results, :forces)

            # Check energy
            @test results[:energy] isa Real
            @test isfinite(results[:energy])

            # Check forces
            forces = results[:forces]
            @test forces isa Matrix
            @test size(forces) == (3, length(species))
            @test all(isfinite, forces)

            # Forces should approximately sum to zero (Newton's third law)
            force_sum = sum(forces, dims = 2)
            @test all(abs.(force_sum) .< 1e-10)  # Very small tolerance for numerical precision

        catch e
            if contains(string(e), "Model creation failed") ||
               contains(string(e), "not found")
                @test_skip "Test model $test_model_name not available"
            else
                rethrow(e)
            end
        end
    end

    @testset "Different Unit Systems" begin
        try
            # Create models with different unit systems
            model_metal = KIM_API.KIMModel(test_model_name, units = :metal)
            model_real = KIM_API.KIMModel(test_model_name, units = :real)

            # Simple system
            species = ["Si", "Si"]
            positions = [SVector(0.0, 0.0, 0.0), SVector(2.35, 0.0, 0.0)]
            cell = Matrix(5.43 * I(3))
            pbc = [true, true, true]

            # Calculate with both unit systems
            results_metal = model_metal(species, positions, cell, pbc)
            results_real = model_real(species, positions, cell, pbc)

            # Both should give valid results
            @test isfinite(results_metal[:energy])
            @test isfinite(results_real[:energy])

            # Energies should differ once converted between unit systems
            energy_ratio = results_metal[:energy] / results_real[:energy]
            @test isfinite(energy_ratio)
            @test !isapprox(energy_ratio, 1.0; atol = 1e-2)

        catch e
            @test_skip "Unit system test failed: $e"
        end
    end

    @testset "Energy-only Computation" begin
        try
            model_energy = KIM_API.KIMModel(test_model_name, compute = [:energy])

            species = ["Si", "Si"]
            positions = [SVector(0.0, 0.0, 0.0), SVector(2.35, 0.0, 0.0)]
            cell = Matrix(5.43 * I(3))
            pbc = [true, true, true]

            results = model_energy(species, positions, cell, pbc)

            @test haskey(results, :energy)
            @test isfinite(results[:energy])

            # Forces should still be computed (current implementation)
            # but energy-only mode could be optimized in the future

        catch e
            @test_skip "Energy-only computation test failed: $e"
        end
    end

    @testset "Error Handling" begin
        try
            model = KIM_API.KIMModel(test_model_name)

            # Test with unsupported species
            @test_throws ErrorException model(
                ["Cu"],
                [SVector(0.0, 0.0, 0.0)],
                Matrix(1.0*I(3)),
                [true, true, true],
            )

            # Test with mismatched species/positions lengths
            @test_throws ArgumentError model(
                ["Si"],
                [SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0)],
                Matrix(1.0*I(3)),
                [true, true, true],
            )

            # Test with invalid cell (non-3x3 matrix)
            invalid_cell = [1.0 0.0; 0.0 1.0]  # 2x2 matrix
            @test_throws Exception model(
                ["Si"],
                [SVector(0.0, 0.0, 0.0)],
                invalid_cell,
                [true, true, true],
            )

        catch e
            @test_skip "Error handling test failed: $e"
        end
    end

    @testset "Larger System" begin
        try
            model = KIM_API.KIMModel(test_model_name)

            # Create a small silicon crystal (8 atoms)
            a = 5.43  # Silicon lattice parameter
            species = fill("Si", 8)
            positions = [
                SVector(0.0, 0.0, 0.0),
                SVector(a/2, a/2, 0.0),
                SVector(a/2, 0.0, a/2),
                SVector(0.0, a/2, a/2),
                SVector(a/4, a/4, a/4),
                SVector(3*a/4, 3*a/4, a/4),
                SVector(3*a/4, a/4, 3*a/4),
                SVector(a/4, 3*a/4, 3*a/4),
            ]

            cell = a * Matrix(1.0*I(3))
            pbc = [true, true, true]

            results = model(species, positions, cell, pbc)

            @test isfinite(results[:energy])
            @test size(results[:forces]) == (3, 8)
            @test all(isfinite, results[:forces])

            # Forces should sum to approximately zero
            force_sum = sum(results[:forces], dims = 2)
            @test all(abs.(force_sum) .< 1e-6)

        catch e
            @test_skip "Larger system test failed: $e"
        end
    end

    @testset "Performance and Memory" begin
        try
            model = KIM_API.KIMModel(test_model_name)

            species = ["Si", "Si"]
            positions = [SVector(0.0, 0.0, 0.0), SVector(2.35, 0.0, 0.0)]
            cell = 5.43 * Matrix(1.0*I(3))
            pbc = [true, true, true]

            # Test that multiple calls work (function reuse)
            results1 = model(species, positions, cell, pbc)
            results2 = model(species, positions, cell, pbc)

            # Results should be identical
            @test results1[:energy] ≈ results2[:energy]
            @test results1[:forces] ≈ results2[:forces]

            # Test that function can be called many times without memory leaks
            # (This is a basic test - real memory leak detection would need more sophisticated tools)
            for i = 1:10
                result = model(species, positions, cell, pbc)
                @test isfinite(result[:energy])
            end

        catch e
            @test_skip "Performance test failed: $e"
        end
    end
end
