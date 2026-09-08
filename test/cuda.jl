# Run explicitly in an environment with CUDA and this checkout available.
using CUDA, EnsembleMCMC, Random, Random123, Test
CUDA.allowscalar(false)
test_rng() = Philox4x((573, 19))
include("kernel.jl")
test_kernel_executor(CuArray)

@testset "CUDA storage and launch boundaries" begin
    initial = randn(MersenneTwister(31), 33, 132)
    scalar(x) = -sum(abs2, x) / 2
    batch!(v, x) = (v .= vec(-sum(abs2, x; dims=1) / 2))
    target = BatchedLogDensity(scalar, batch!)
    state = initialize(test_rng(), target, CuArray(initial); move=DEMove(), executor=KernelExecutor())
    expected = sample!(initialize(test_rng(), scalar, initial; move=DEMove()), 2)
    draws = sample!(state, 2)
    @test Array(draws.positions) ≈ expected.positions
    @test Array(draws.accepted) == expected.accepted
    saved = snapshot(state)
    @test all(x -> x isa CUDA.AnyCuArray,
        (draws.positions, draws.logdensities, draws.accepted, saved.positions[1], saved.logdensities, saved.accepted))

    event = CUDA.CuEvent()
    throws!(v, x) = (fill!(v, 0); CUDA.record(event); error("callback failed"))
    failed = initialize(test_rng(), BatchedLogDensity(scalar, throws!), CuArray(initial);
        executor=KernelExecutor())
    @test_throws ErrorException step!(failed)
    @test CUDA.isdone(event)
end
