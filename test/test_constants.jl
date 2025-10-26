# test_constants.jl - Test KIM_API constants and enumerations

@testset "Constants and Enumerations" begin

    @testset "Enumeration Values" begin
        # Test that enums have correct integer values (matching KIM-API C constants)
        @test Int(KIM_API.zeroBased) == 0
        @test Int(KIM_API.oneBased) == 1

        @test Int(KIM_API.A) == 1
        @test Int(KIM_API.Bohr) == 2
        @test Int(KIM_API.eV) == 3
        @test Int(KIM_API.e) == 2
        @test Int(KIM_API.K) == 1
        @test Int(KIM_API.ps) == 2

        @test Int(KIM_API.numberOfParticles) == 0
        @test Int(KIM_API.particleSpeciesCodes) == 1
        @test Int(KIM_API.coordinates) == 3
        @test Int(KIM_API.partialEnergy) == 4
        @test Int(KIM_API.partialForces) == 5

        @test Int(KIM_API.GetNeighborList) == 0
        @test Int(KIM_API.c) == 1
        @test Int(KIM_API.required) == 2
        @test Int(KIM_API.optional) == 3
        @test Int(KIM_API.notSupported) == 1
    end

    @testset "Unit Styles" begin
        # Test predefined unit styles
        metal = KIM_API.UNIT_STYLES.metal
        @test metal.length == KIM_API.A
        @test metal.energy == KIM_API.eV
        @test metal.charge == KIM_API.e
        @test metal.temperature == KIM_API.K
        @test metal.time == KIM_API.ps

        real = KIM_API.UNIT_STYLES.real
        @test real.length == KIM_API.A
        @test real.energy == KIM_API.kcal_mol
        @test real.time == KIM_API.fs

        si = KIM_API.UNIT_STYLES.si
        @test si.length == KIM_API.m
        @test si.energy == KIM_API.J
        @test si.charge == KIM_API.C
        @test si.time == KIM_API.s

        electron = KIM_API.UNIT_STYLES.electron
        @test electron.length == KIM_API.Bohr
        @test electron.energy == KIM_API.Hartree
        @test electron.time == KIM_API.fs
    end

    @testset "Unit Style Function" begin
        # Test get_lammps_style_units function
        metal_units = KIM_API.get_lammps_style_units(:metal)
        @test metal_units.length == KIM_API.A
        @test metal_units.energy == KIM_API.eV

        real_units = KIM_API.get_lammps_style_units(:real)
        @test real_units.energy == KIM_API.kcal_mol

        si_units = KIM_API.get_lammps_style_units(:si)
        @test si_units.length == KIM_API.m
        @test si_units.energy == KIM_API.J

        # Test error for invalid style
        @test_throws ErrorException KIM_API.get_lammps_style_units(:invalid)
    end

    if isdefined(KIM_API, :libkim) && KIM_API.libkim != ""
        @testset "KIM-API Constant Functions" begin
            # Test constant lookup functions (require KIM-API library)
            @test KIM_API.get_numbering("zeroBased") == 0
            @test KIM_API.get_numbering("oneBased") == 1

            @test KIM_API.get_length_unit("A") == 1
            @test KIM_API.get_length_unit("Bohr") == 2

            @test KIM_API.get_energy_unit("eV") == 3
            @test KIM_API.get_energy_unit("J") == 5

            @test KIM_API.get_charge_unit("e") == 2
            @test KIM_API.get_charge_unit("C") == 1

            @test KIM_API.get_temperature_unit("K") == 1

            @test KIM_API.get_time_unit("fs") == 1
            @test KIM_API.get_time_unit("ps") == 2
        end

        @testset "Reverse Lookup Functions" begin
            # Test constant to string conversion
            @test KIM_API.numbering_to_string(0) == "zeroBased"
            @test KIM_API.numbering_to_string(1) == "oneBased"

            @test KIM_API.length_unit_to_string(1) == "A"
            @test KIM_API.length_unit_to_string(2) == "Bohr"

            @test KIM_API.energy_unit_to_string(3) == "eV"
            @test KIM_API.energy_unit_to_string(5) == "J"

            # Test unknown values
            @test KIM_API.numbering_to_string(999) == "unknown"
            @test KIM_API.length_unit_to_string(-1) == "unknown"
        end
    end
end
