@testset "Ham-terms" begin
    @testset "Constructor" begin
        HoppingDiagonal(1.0, 0, [0, 0])
        HoppingDiagonal{Float64}(1.0, 0, [0, 0])

        HoppingOffdiagonal(1.0, 0, 1, [0, 0], [0, 0])
        HoppingOffdiagonal{Float64}(1.0, 0, 1, [0, 0], [0, 0])
        HoppingOffdiagonal{ComplexF64}(1.0+0.0im, 0, 1, [0, 0], [0, 0])

        HoppingOffdiagonal(1.0, 0, 0, [0, 0], [1, 0])
        @test_throws ArgumentError HoppingOffdiagonal(1.0, 0, 0, [0, 0], [0, 0])

    end

    @testset "Constructor: Exceptions" begin
        @test_throws ArgumentError HoppingOffdiagonal(1.0, 0, 0, [0,], [0, 0])          # Dimension matching
        @test_throws ArgumentError HoppingOffdiagonal(1.0, 0, 0, [0, 0], [0, 0])        # Not Offdiagonal
        @test_throws ArgumentError InteractionDiagonal(1.0, 0, 0, [0,], [0, 0])         # Dimension matching
        @test_throws ArgumentError InteractionDiagonal(1.0, 0, 0, [0, 0], [0, 0])       # Fermi Statistics

        #@test_throws ArgumentError InteractionOffdiagonal(1.0, 0, 0, 0, 0, [0,], [0, 0])
        # Fermi Statistics
        @test_throws ArgumentError InteractionOffdiagonal(1.0, 0, 0, 1, 2, [0, 0], [0, 0], [0, 0], [0, 0])
        @test_throws ArgumentError InteractionOffdiagonal(1.0, 1, 2, 0, 0, [0, 0], [0, 0], [0, 0], [0, 0])

        # Not Offdiagonal
        @test_throws ArgumentError InteractionOffdiagonal(1.0, 0, 1, 0, 1, [0, 0], [0, 0], [0, 0], [0, 0])
        InteractionOffdiagonal(1.0, 0, 1, 0, 2, [0, 0], [0, 0], [0, 0], [0, 0])        # Orbital difference
        InteractionOffdiagonal(1.0, 0, 1, 0, 1, [0, 0], [0, 0], [0, 0], [1, 0])        # Unitcell difference
    end

    @testset "Ordering" begin
        let
            ho1 = HoppingOffdiagonal(1.0 + 0.5im, 0, 1, [0, 0], [1, 0])
            ho2 = HoppingOffdiagonal(1.0 + 0.5im, 1, 0, [1, 0], [0, 0])
            ho3 = HoppingOffdiagonal(1.0 - 0.5im, 1, 0, [1, 0], [0, 0])

            @test !isapprox(ho1, ho2)
            @test isapprox(ho1, ho3)
        end

        let
            io1 = InteractionDiagonal(1.0, 0, 1, [0, 0], [1, 0])
            io2 = InteractionDiagonal(1.0, 1, 0, [1, 0], [0, 0])

            @test isapprox(io1, io2)
        end
        # TODO(kyungminlee, 20190102): Keep writing the testset
    end

    @testset "local" begin
    end

    @testset "isapprox" begin
    end
end
