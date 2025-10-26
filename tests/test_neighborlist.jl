# AI GENERATED TESTS
# TODO: Review and clean

using Test
using Random
using StaticArrays
using KIM_API

# Include the neighbors.jl file if it's not part of module
# include("../src/neighbors.jl")

@testset "KIM Neighbors Tests" begin

    @testset "NeighborListContainer" begin
        # Test basic construction
        dummy_fn = x -> [1, 2, 3]
        container = KIM_API.NeighborListContainer(dummy_fn)

        @test container.get_neighbors === dummy_fn
        @test isempty(container.neighbor_storage)
        @test eltype(container.neighbor_storage) == Int32
    end

    @testset "create_neighborlist_closure" begin
        # Simple 2D system for easier testing
        positions = [
            SVector(0.0, 0.0, 0.0),
            SVector(1.0, 0.0, 0.0),
            SVector(0.0, 1.0, 0.0),
            SVector(1.0, 1.0, 0.0),
        ]

        # Test without PBC
        cell = [10.0 0.0 0.0; 0.0 10.0 0.0; 0.0 0.0 10.0]
        pbc = [false, false, false]
        cutoff = 1.5

        get_neighbors = KIM_API.create_neighborlist_closure(positions, cutoff, cell, pbc)

        # Particle 1 should have neighbors 2 and 3 (distance = 1.0)
        neighbors = get_neighbors(1)
        @test length(neighbors) == 2
        @test 2 in neighbors
        @test 3 in neighbors

        # Test error when cell is missing with PBC
        @test_throws ArgumentError KIM_API.create_neighborlist_closure(
            positions,
            cutoff,
            nothing,
            [true, true, true],
        )
    end

    @testset "create_kim_neighborlist_dataobject" begin
        positions = [SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0)]
        cutoffs = [1.5, 2.5, 3.5]
        cell = [10.0 0.0 0.0; 0.0 10.0 0.0; 0.0 0.0 10.0]

        nl_containers = KIM_API.create_kim_neighborlist_dataobject(positions, cutoffs, cell)

        @test length(nl_containers) == 3
        @test all(c isa KIM_API.NeighborListContainer for c in nl_containers)
    end

    @testset "kim_neighbors_function" begin
        # Create test data
        positions = [SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(2.0, 0.0, 0.0)]
        cutoffs = [1.5, 2.5]
        cell = [10.0 0.0 0.0; 0.0 10.0 0.0; 0.0 0.0 10.0]

        nl_containers = KIM_API.create_kim_neighborlist_dataobject(
            positions,
            cutoffs,
            cell,
            [false, false, false],
        )
        dataObject = pointer_from_objref(nl_containers)

        # Test valid neighbor list request
        numberOfNeighbors = Ref{Cint}(0)
        neighborsPtr = Ref{Ptr{Cint}}(C_NULL)

        # Request neighbors for particle 0 (0-based) with cutoff index 0
        result = KIM_API.kim_neighbors_function(
            dataObject,
            Cint(2),      # numberOfNeighborLists
            C_NULL,       # cutoffs (not used in current implementation)
            Cint(0),      # neighborListIndex
            Cint(0),      # particleNumber (0-based)
            numberOfNeighbors,
            neighborsPtr,
        )

        @test result == Cint(0)  # Success
        @test numberOfNeighbors[] == 1  # Should have 1 neighbor

        # Check the actual neighbors
        neighbors = unsafe_wrap(Array, neighborsPtr[], numberOfNeighbors[])
        @test neighbors[1] == 1  # Particle 1 (0-based)

        # Test with larger cutoff (index 1)
        result = KIM_API.kim_neighbors_function(
            dataObject,
            Cint(2),
            C_NULL,
            Cint(1),      # neighborListIndex = 1
            Cint(0),      # particleNumber = 0
            numberOfNeighbors,
            neighborsPtr,
        )

        @test result == Cint(0)
        @test numberOfNeighbors[] == 2  # Should have 2 neighbors

        # Test invalid neighbor list index
        result = KIM_API.kim_neighbors_function(
            dataObject,
            Cint(2),
            C_NULL,
            Cint(2),      # Invalid index
            Cint(0),
            numberOfNeighbors,
            neighborsPtr,
        )

        @test result == Cint(1)  # Error

        # Test negative index
        result = KIM_API.kim_neighbors_function(
            dataObject,
            Cint(2),
            C_NULL,
            Cint(-1),     # Invalid index
            Cint(0),
            numberOfNeighbors,
            neighborsPtr,
        )

        @test result == Cint(1)  # Error
    end

    @testset "Index conversion" begin
        # Test that KIM's 0-based indices are properly converted
        positions = [SVector(float(i), 0.0, 0.0) for i = 0:4]
        cutoff = 1.5
        cell = [10.0 0.0 0.0; 0.0 10.0 0.0; 0.0 0.0 10.0]

        get_neighbors =
            KIM_API.create_neighborlist_closure(positions, cutoff, cell, [false, false, false])

        # Julia: particle 2 (1-based) should have neighbors 1 and 3
        julia_neighbors = get_neighbors(2)
        @test 1 in julia_neighbors
        @test 3 in julia_neighbors

        # Through KIM interface (0-based)
        nl_containers = [KIM_API.NeighborListContainer(get_neighbors)]
        dataObject = pointer_from_objref(nl_containers)

        numberOfNeighbors = Ref{Cint}(0)
        neighborsPtr = Ref{Ptr{Cint}}(C_NULL)

        result = KIM_API.kim_neighbors_function(
            dataObject,
            Cint(1),
            C_NULL,
            Cint(0),
            Cint(1),      # Particle 1 in 0-based = particle 2 in 1-based
            numberOfNeighbors,
            neighborsPtr,
        )

        @test result == Cint(0)
        kim_neighbors = unsafe_wrap(Array, neighborsPtr[], numberOfNeighbors[])
        @test 0 in kim_neighbors  # 0-based index for particle 1
        @test 2 in kim_neighbors  # 0-based index for particle 3
    end

    @testset "Large system performance" begin
        # Test with larger system to ensure it scales
        Random.seed!(42)
        n_atoms = 1000
        positions = [SVector(10*rand(), 10*rand(), 10*rand()) for _ = 1:n_atoms]
        cutoff = 1.0
        cell = [10.0 0.0 0.0; 0.0 10.0 0.0; 0.0 0.0 10.0]

        nl_containers = KIM_API.create_kim_neighborlist_dataobject(
            positions,
            [cutoff],
            cell,
            [true, true, true],
        )

        @test length(nl_containers) == 1

        # Just verify it runs without error
        get_neighbors = nl_containers[1].get_neighbors
        neighbors = get_neighbors(1)
        @test neighbors isa Vector
    end
end
