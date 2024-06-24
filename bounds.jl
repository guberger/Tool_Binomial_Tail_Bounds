using Distributions

"""
    binomial_tail(N, p, ϵ; method=:exact)

Bound on
```math
P[\\frac{1}{N} \\sum_{i=1}^N X_i - p > \\epsilon]
```
where ``X_i`` is a Bernoulli variable with ``P[X_i=1] = p``.
"""
function binomial_tail(N, p, ϵ)
    nothing
end

function check_args(N, p, ϵ)
    @assert isinteger(N) && N > 0
    @assert zero(p) ≤ p ≤ one(p)
    @assert zero(ϵ) ≤ ϵ ≤ one(ϵ)
    nothing
end

function exact_value(N, p, ϵ)
    check_args(N, p, ϵ)
    D = Binomial(N, p)
    k = floor(Int, (p + ϵ) * N + 1)
    return float(ccdf(D, k))
end

function hoeffding_bound(N, ::Any, ϵ)
    check_args(N, 1, ϵ)
    return exp(-2 * N * ϵ^2)
end

function bernstein_bound(N, p, ϵ)
    check_args(N, p, ϵ)
    if iszero(p) || isone(p)
        return 0.0
    end
    # here 0 < p < 1
    σ2 = p * (1 - p) # σ2 > 0
    D = 2 * (σ2 + ϵ / 3)
    return exp(-N * ϵ^2 / D)
end

function bennett_bound(N, p, ϵ)
    check_args(N, p, ϵ)
    if iszero(p) || isone(p)
        return 0.0
    end
    # here 0 < p < 1
    σ2 = p * (1 - p) # σ2 > 0
    r = ϵ / σ2
    θ = (1 + r) * log(1 + r) - r
    return exp(-N * σ2 * θ)
end

function sanov_bound(N, p, ϵ)
    check_args(N, p, ϵ)
    if iszero(p) || isone(p) || p + ϵ ≥ 1
        return 0.0
    end
    # here 0 < p ≤ p + ϵ < 1
    q = p + ϵ
    D = q * log(q / p) + (1 - q) * log((1 - q) / (1 - p))
    return exp(-N * D)
end