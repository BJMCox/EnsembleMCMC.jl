abstract type AbstractExecutor end
"""
    SequentialExec()

Evaluate walker proposals serially within each group of a sweep.
"""
struct SequentialExec <: AbstractExecutor end
"""
    MultiThreadedExec()

Evaluate walker proposals with Julia threads within each group of a sweep.
Groups still advance in order. The log-density function must support concurrent
calls. A deterministic target gives the same trajectory as [`SequentialExec`](@ref)
for the same initial state, move, RNG, and walker IDs.
"""
struct MultiThreadedExec <: AbstractExecutor end

"""
    MoveMixture(moves, weights)

Select one move per complete sweep. Integer weights give a repeating weighted
cycle in component order. Other real weights give fixed random selection
probabilities after normalization. Weights must be finite and nonnegative,
with at least one positive weight.
"""
struct MoveMixture{M<:Tuple,W<:AbstractVector}
    moves::M
    weights::W
    function MoveMixture(moves, weights::AbstractVector{<:Real})
        Base.require_one_based_indexing(weights)
        ms = Tuple(moves)
        !isempty(ms) && all(m -> m isa AbstractEnsembleMove, ms) ||
            throw(ArgumentError("A mixture requires ensemble moves"))
        length(ms) == length(weights) || throw(DimensionMismatch("Moves and weights must match"))
        all(w -> isfinite(w) && w >= 0, weights) && any(>(0), weights) ||
            throw(ArgumentError("Weights must be finite, nonnegative and have positive mass"))
        ws = if eltype(weights) <: Integer
            total = sum(big, weights)
            total <= typemax(Int) || throw(ArgumentError("Mixture cycle is too long"))
            Int.(weights)
        else
            values = collect(promote(map(float, weights)...))
            scaled = values ./ maximum(values)
            scaled ./ sum(scaled)
        end
        new{typeof(ms),typeof(ws)}(ms, ws)
    end
end

const _PROPOSALS_PER_PURPOSE = typemax(Int16) - 2
const _COMPANION_PURPOSE = 5
const _SCALE_PURPOSE = 6

_stream_index(purpose, proposal) = (purpose - 1) * _PROPOSALS_PER_PURPOSE + proposal
_walker_rngpart(part, purpose, proposal) = RNGPartition(
    AbstractRNG(part, _stream_index(purpose, proposal)), Base.OneTo(typemax(Int32) - 2),
)

mutable struct EnsembleState{F,M,W,E,R,P,V,L}
    logdensity::F
    moves::M
    weights::W
    executor::E
    rng::R
    cycle_partition::P
    walker_rngs::Vector{R}
    walker_ids::Vector{Int}
    walker_order::Vector{Int}
    positions::V
    candidates::V
    logdensities::L
    candidate_logdensities::L
    accepted::Vector{Bool}
    attempts::Vector{Int}
    accepts::Vector{Int}
    step::Int
    active_index::Int
    valid::Bool
end

"""
    initialize(rng, logdensity, initial; move=StretchMove(), executor=SequentialExec(),
               walker_ids=eachindex(initial))

Create one ensemble from a vector of finite coordinate vectors. The coordinates
must have full affine rank and finite initial log densities. The state owns
copies of coordinates and RNG. Supported RNGs are Random123 Philox4x/Threefry4x
with UInt64 counters.
Use independent RNG seeds or streams for independent ensembles.

Walker IDs remain fixed throughout the run. The density must not mutate its input,
and must support concurrent calls with MultiThreadedExec. Do not mutate state fields.
"""
function initialize(
    rng::Union{Philox4x{UInt64},Threefry4x{UInt64}}, logdensity, initial::AbstractVector;
    move = StretchMove(), executor::AbstractExecutor = SequentialExec(),
    walker_ids = collect(eachindex(initial)),
)
    isempty(initial) && throw(ArgumentError("An ensemble cannot be empty"))
    Base.require_one_based_indexing(initial)
    d = length(first(initial))
    d > 0 && all(x -> length(x) == d, initial) ||
        throw(DimensionMismatch("Walker dimensions must agree and be positive"))
    T = float(promote_type(map(eltype, initial)...))
    T <: AbstractFloat || throw(ArgumentError("Coordinates must be real floating-point values"))
    positions = [Vector{T}(x) for x in initial]
    all(x -> all(isfinite, x), positions) || throw(ArgumentError("Coordinates must be finite"))
    moves = move isa MoveMixture ? map(m -> prepare_move(m, T, d), move.moves) :
        (prepare_move(move, T, d),)
    length(moves) <= _PROPOSALS_PER_PURPOSE || throw(ArgumentError("Too many mixture components"))
    n = length(positions)
    n >= maximum(m -> minimum_walkers(m, d), moves) ||
        throw(ArgumentError("Too few walkers for the dimension and selected moves"))
    centered = reduce(hcat, [x .- first(positions) for x in positions])
    rank(centered; rtol = max(size(centered)...) * eps(T)) == d ||
        throw(ArgumentError("Initial walkers must have full affine rank"))
    ids = collect(Int, walker_ids)
    length(ids) == n && length(unique(ids)) == n &&
        all(id -> 1 <= id <= typemax(Int32) - 2, ids) ||
        throw(ArgumentError("Walker IDs must be distinct positive stream indices"))
    logds = collect(promote(map(x -> _checked_logdensity(logdensity, x), positions)...))
    all(isfinite, logds) || throw(ArgumentError("Initial log densities must be finite"))
    owned_rng = copy(rng)
    cycle_part = RNGPartition(owned_rng, 0:(typemax(Int16) - 2))
    weights = move isa MoveMixture ? copy(move.weights) : nothing
    return EnsembleState(
        logdensity, moves, weights, executor, owned_rng, cycle_part,
        [rngpart_createrng(typeof(rng)) for _ in 1:n], ids, sortperm(ids),
        positions, deepcopy(positions), logds, copy(logds), fill(false, n),
        zeros(Int, length(moves)), zeros(Int, length(moves)), 0, 1, true,
    )
end

function _checked_logdensity(f, x)
    value = f(x)
    value isa Real || throw(ArgumentError("Log density must return a real number"))
    (isnan(value) || value == Inf) && throw(DomainError(value, "Invalid log density"))
    return float(value)
end

_select_move(::Nothing, rng, step) = 1
function _select_move(weights::Vector{Int}, rng, step)
    offset = mod1(step, sum(weights))
    total = 0
    for i in eachindex(weights)
        total += weights[i]
        offset <= total && return i
    end
    error("Invalid mixture cycle")
end
function _select_move(weights::AbstractVector{T}, rng, step) where {T<:AbstractFloat}
    u = rand(rng, T) * sum(weights)
    total = zero(T)
    for i in eachindex(weights)
        total += weights[i]
        u < total && return i
    end
    return findlast(>(0), weights)
end

function _groups(move, part, proposal_idx, order)
    rng = AbstractRNG(part, _stream_index(4, proposal_idx))
    permutation = randperm(rng, length(order))
    ngroups = group_count(move)
    if ngroups == 2
        split = fld(length(order), 2)
        left = order[permutation[begin:split]]
        right = order[permutation[(split + 1):end]]
        return rand(rng, Bool) ? (left, right) : (right, left)
    end
    size, extra = divrem(length(order), ngroups)
    groups = Vector{Vector{Int}}(undef, ngroups)
    start = 1
    for i in eachindex(groups)
        stop = start + size - 1 + (i <= extra)
        groups[i] = order[permutation[start:stop]]
        start = stop + 1
    end
    shuffle!(rng, groups)
    return groups
end

_complement(groups::Tuple, i) = i == 1 ? (groups[2],) : (groups[1],)
_complement(groups::AbstractVector, i) = groups[eachindex(groups) .!= i]

function _evaluate!(state, move, part, proposal_idx, i, complement, acceptance_part)
    rng = state.walker_rngs[i]
    log_hastings = propose!(state.candidates[i], move, state.positions, i,
        complement, rng, part, proposal_idx, state.walker_ids[i])
    if isnothing(log_hastings)
        state.accepted[i] = false
        return nothing
    end
    logd = convert(eltype(state.logdensities), _checked_logdensity(state.logdensity, state.candidates[i]))
    (isnan(logd) || logd == Inf) && throw(DomainError(logd, "Invalid stored log density"))
    state.candidate_logdensities[i] = logd
    T = eltype(state.positions[i])
    logratio = convert(T, log_hastings + logd - state.logdensities[i])
    probability = isnan(logratio) ? zero(T) : clamp(exp(logratio), zero(T), one(T))
    set_rng!(rng, acceptance_part, state.walker_ids[i])
    state.accepted[i] = rand(rng, T) < probability
    return nothing
end

function _evaluate_group!(::SequentialExec, state, move, part, idx, group, complement, acceptance_part)
    for i in group
        _evaluate!(state, move, part, idx, i, complement, acceptance_part)
    end
end
function _evaluate_group!(::MultiThreadedExec, state, move, part, idx, group, complement, acceptance_part)
    Threads.@threads for k in eachindex(group)
        _evaluate!(state, move, part, idx, group[k], complement, acceptance_part)
    end
end

"""
    step!(state)

Advance one full ensemble sweep. Mutate and return the state.
An exception during a sweep invalidates the state. Initialize a fresh state
after correcting the target instead of resuming a partially completed sweep.
"""
function step!(state::EnsembleState)
    state.valid || throw(ArgumentError("Cannot resume a state after a failed sweep"))
    state.step < typemax(Int32) - 2 || throw(ArgumentError("RNG step range exhausted"))
    state.valid = false
    set_rng!(state.rng, state.cycle_partition, 0)
    steps = RNGPartition(state.rng, 0:(typemax(Int32) - 2))
    set_rng!(state.rng, steps, state.step)
    part = RNGPartition(state.rng, Base.OneTo(6 * _PROPOSALS_PER_PURPOSE))
    selection_rng = AbstractRNG(part, _stream_index(1, 1))
    idx = _select_move(state.weights, selection_rng, state.step + 1)
    move = state.moves[idx]
    groups = _groups(move, part, idx, state.walker_order)
    acceptance_part = _walker_rngpart(part, 3, idx)
    for active in eachindex(groups)
        group = groups[active]
        _evaluate_group!(state.executor, state, move, part, idx, group,
            _complement(groups, active), acceptance_part)
        for i in group
            if state.accepted[i]
                copyto!(state.positions[i], state.candidates[i])
                state.logdensities[i] = state.candidate_logdensities[i]
            end
        end
    end
    state.active_index = idx
    state.attempts[idx] += length(state.positions)
    state.accepts[idx] += count(state.accepted)
    state.step += 1
    state.valid = true
    return state
end

"""
    sample!(state, nsweeps)

Advance and collect complete sweeps. Positions have axes (coordinate, walker,
sweep). Rejections repeat the current state. Returned arrays own their storage.
Each state represents one coupled ensemble, not independent walker chains.

Return a named tuple with `positions`, `logdensities` and `accepted` (both indexed
by walker and sweep), `proposal_indices` (one move index per sweep), and
`walker_ids` (one ID per walker). The initial positions are not included.
"""
function sample!(state::EnsembleState, nsweeps::Integer)
    nsweeps >= 0 || throw(ArgumentError("Sweep count must be nonnegative"))
    n = length(state.positions)
    d = length(first(state.positions))
    positions = Array{eltype(first(state.positions))}(undef, d, n, nsweeps)
    logdensities = Matrix{eltype(state.logdensities)}(undef, n, nsweeps)
    accepted = Matrix{Bool}(undef, n, nsweeps)
    proposal_indices = Vector{Int}(undef, nsweeps)
    for sweep in 1:nsweeps
        step!(state)
        for i in 1:n
            positions[:, i, sweep] .= state.positions[i]
        end
        logdensities[:, sweep] .= state.logdensities
        accepted[:, sweep] .= state.accepted
        proposal_indices[sweep] = state.active_index
    end
    return (; positions, logdensities, accepted, proposal_indices,
        walker_ids = copy(state.walker_ids))
end
