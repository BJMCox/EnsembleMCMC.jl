using EnsembleMCMC
using LinearAlgebra
using Random
using Random123
using Statistics
using Test

gaussian_logdensity(x) = -sum(abs2, x) / 2
initial_walkers() = [randn(Philox4x((731, i)), 2) for i in 1:24]
test_rng() = Philox4x((573, 19))

@testset "EnsembleMCMC" begin
    @testset "Gaussian target: $(typeof(move))" for move in
        (StretchMove(), DEMove(), DESnookerMove())
        state = initialize(test_rng(), gaussian_logdensity, initial_walkers(); move)
        sample!(state, 500)
        draws = sample!(state, 2_000)
        coordinates = reshape(draws.positions, 2, :)
        @test maximum(abs, vec(mean(coordinates; dims=2))) < 0.12
        @test cov(coordinates; dims=2) ≈ Matrix{Float64}(I, 2, 2) atol=0.15 rtol=0
    end

    @testset "Affine trajectory: $(typeof(move))" for move in
        (StretchMove(), DEMove())
        matrix = [1.3 0.4; -0.2 0.8]
        shift = [-0.6, 1.1]
        initial = initial_walkers()
        target(y) = gaussian_logdensity(matrix \ (y - shift))
        reference = initialize(test_rng(), gaussian_logdensity, initial; move)
        transformed = initialize(test_rng(), target, [matrix * x + shift for x in initial]; move)
        expected = sample!(reference, 32)
        actual = sample!(transformed, 32)
        @test actual.accepted == expected.accepted
        @test all(
            actual.positions[:, w, s] ≈ matrix * expected.positions[:, w, s] + shift
            for w in axes(expected.positions, 2), s in axes(expected.positions, 3)
        )
    end

    @testset "Thread and walker order replay" begin
        initial = initial_walkers()
        permutation = reverse(eachindex(initial))
        move = MoveMixture((StretchMove(), DEMove(), DESnookerMove()), [1, 1, 1])
        reference = initialize(test_rng(), gaussian_logdensity, initial; move)
        threaded = initialize(test_rng(), gaussian_logdensity, initial;
            move, executor=MultiThreadedExec())
        reordered = initialize(test_rng(), gaussian_logdensity, initial[permutation];
            move, walker_ids=collect(permutation))
        expected = sample!(reference, 30)
        actual = sample!(threaded, 30)
        reordered_draws = sample!(reordered, 30)
        @test actual == expected
        @test reordered_draws.positions[:, permutation, :] == expected.positions
        @test reordered_draws.logdensities[permutation, :] == expected.logdensities
        @test reordered_draws.accepted[permutation, :] == expected.accepted
        @test reordered_draws.walker_ids[permutation] == expected.walker_ids
    end

    @testset "Resume and snapshot ownership" begin
        initial = initial_walkers()
        saved_initial = deepcopy(initial)
        move = MoveMixture((StretchMove(), DEMove()), [2, 1])
        state = initialize(test_rng(), gaussian_logdensity, initial; move)
        reference = initialize(test_rng(), gaussian_logdensity, saved_initial; move)
        initial[1][1] = 1e6
        first_part = sample!(state, 7)
        second_part = sample!(state, 11)
        expected = sample!(reference, 18)
        @test cat(first_part.positions, second_part.positions; dims=3) == expected.positions
        @test vcat(first_part.proposal_indices, second_part.proposal_indices) == expected.proposal_indices
        @test hcat(first_part.logdensities, second_part.logdensities) == expected.logdensities
        @test hcat(first_part.accepted, second_part.accepted) == expected.accepted
        snapshot = deepcopy(second_part)
        @test step!(state) === state
        @test second_part == snapshot
        second_part.positions[1, 1, end] = 1e6
        @test state.positions == step!(reference).positions
        @test state.step == 19
    end

    @testset "Fixed mixture selection" begin
        moves = (StretchMove(), DEMove())
        deterministic = initialize(test_rng(), gaussian_logdensity, initial_walkers();
            move=MoveMixture(moves, [2, 1]))
        draws = sample!(deterministic, 9)
        @test draws.proposal_indices == repeat([1, 1, 2], 3)
        @test deterministic.attempts == [6, 3] .* length(initial_walkers())
        @test sum(deterministic.accepts) == count(draws.accepted)
        random = initialize(test_rng(), gaussian_logdensity, initial_walkers();
            move=MoveMixture(moves, [0.7, 0.3]))
        random_draws = sample!(random, 100)
        @test 50 < count(==(1), random_draws.proposal_indices) < 90
        @test random.attempts ==
            [count(==(i), random_draws.proposal_indices) for i in 1:2] .* length(initial_walkers())
    end

    @testset "Cached densities and degenerate snooker geometry" begin
        calls = Ref(0)
        target(x) = (calls[] += 1; gaussian_logdensity(x))
        initial = initial_walkers()
        state = initialize(test_rng(), target, initial)
        @test calls[] == length(initial)
        step!(state)
        @test calls[] == 2length(initial)
        calls[] = 0
        coincident = [[0.0], [0.0], [0.0], [1.0]]
        degenerate = initialize(test_rng(), target, coincident; move=DESnookerMove())
        step!(degenerate)
        @test calls[] < 2length(coincident)
        @test count(!, degenerate.accepted) > 0
    end

    @testset "Failed sweeps cannot resume" begin
        calls = Ref(0)
        enabled = Ref(false)
        initial = initial_walkers()
        target(x) = begin
            calls[] += 1
            enabled[] && calls[] > length(initial) + length(initial) ÷ 2 && error("target failed")
            gaussian_logdensity(x)
        end
        state = initialize(test_rng(), target, initial)
        enabled[] = true
        @test_throws ErrorException step!(state)
        before_retry = calls[]
        enabled[] = false
        @test_throws ArgumentError step!(state)
        @test calls[] == before_retry
    end
end
