import RigidBodyDynamics: hat, rotation_vector_rate

@testset "util" begin
    @testset "rotation vector rate" begin
        for ϕ in (rand(SVector{3}), zeros(SVector{3})) # exponential coordinates (rotation vector)
            ω = rand(SVector{3}) # angular velocity in body frame
            R = RotMatrix(RodriguesVec(ϕ...))
            Ṙ = R * hat(ω)
            ϕ̇ = rotation_vector_rate(ϕ, ω)
            Θ = norm(ϕ)
            if Θ > eps(Θ)
                ϕ_autodiff = SVector{3}(create_autodiff(ϕ, ϕ̇))
                R_autodiff = RotMatrix(RodriguesVec(ϕ_autodiff...))
                Ṙ_from_autodiff = map(x -> ForwardDiff.partials(x)[1], R_autodiff)
                @test isapprox(Ṙ_from_autodiff, Ṙ)
            else
                @test isapprox(ϕ̇, ω) # limit case; hard to test using autodiff because of division by zero
            end
        end
    end

    @testset "AngleAxis from two vectors" begin
        angle_axis_test(from, to, rot, atol) = isapprox(rot * from * norm(to) / norm(from), to; atol = atol)

        for i = 1 : 10000
            from = rand(SVector{3, Float64})
            to = rand(SVector{3, Float64})
            rot = AngleAxis(from, to)
            @test angle_axis_test(from, to, rot, 1e-10)
        end

        # degenerate cases
        for i = 1 : 10000
            from = rand(SVector{3, Float64})
            to = randn() * from # either from and to are aligned, or in opposite directions
            rot = AngleAxis(from, to)
            @test angle_axis_test(from, to, rot, 1e-7)
        end
        for direction = 1 : 3
            for i = 1 : 10000
                from = @SVector [ifelse(i == direction, 1., 0.) for i = 1 : 3] # unit vector in direction 'direction'
                to = randn() * from
                rot = AngleAxis(from, to)
                @test angle_axis_test(from, to, rot, 1e-7)
            end
        end
        @test_throws ArgumentError AngleAxis(zero(SVector{3}), rand(SVector{3}))
        @test_throws ArgumentError AngleAxis(rand(SVector{3}), zero(SVector{3}))
    end
end
