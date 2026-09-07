@testset "Snooker geometry across coordinate scales" begin
    for T in (Float32, Float64)
        initial = T[-1 0 1 0 -1 -1 1 1; 0 -1 0 1 -1 1 -1 1]
        reference = sample!(initialize(test_rng(), _ -> zero(T), initial;
            move=DESnookerMove()), 1)
        exponent = T === Float32 ? 25 : 200
        for scale in (T(10)^(-exponent), T(10)^exponent)
            draws = sample!(initialize(test_rng(), _ -> zero(T), initial .* scale;
                move=DESnookerMove()), 1)
            @test draws.accepted == reference.accepted
            @test draws.positions ./ scale ≈ reference.positions rtol=100eps(T)
        end
    end
end
