###
# Calculate profiles using the naive, polynomial-time algorithm.
# Used for comparison with more efficient methods.
###
using Combinatorics
using BenchmarkTools
using Base.Threads  

function standardize(v::AbstractVector{Int})
    sorted_vals = sort(v)
    rank = Dict(val => i for (i,val) in enumerate(sorted_vals))
    return [rank[x] for x in v]
end

# cache of permutation index dictionaries by k
const _perm_index_cache = Dict{Int, Dict{Tuple{Vararg{Int}}, Int}}()

# build (or fetch) a mapping from permutation pattern (Tuple) -> lexicographic index
function _get_perm_index(k::Int)
    if haskey(_perm_index_cache, k)
        return _perm_index_cache[k]
    end
    perm_lex = sort(collect(permutations(1:k)))
    d = Dict{Tuple{Vararg{Int}}, Int}()
    for (i,p) in enumerate(perm_lex)
        d[Tuple(p)] = i
    end
    _perm_index_cache[k] = d
    return d
end

"""
    profile_bruteforce(perm::AbstractVector{Int}, k::Int)

Compute the profile vector for subsequence order `k` (2 ≤ k ≤ 6).
Returns a vector of length factorial(k) where entry i counts how many length-k
subsequences of `perm` follow the i-th lexicographic pattern.
"""
function profile_bruteforce(perm::AbstractVector{Int}, k::Int)
    if k < 2 || k > 6
        throw(ArgumentError("k must be between 2 and 6"))
    end
    n = length(perm)
    if n < k
        throw(ArgumentError("input vector length must be at least k"))
    end

    perm_index = _get_perm_index(k)
    freq = zeros(Int, factorial(k))

    for inds in combinations(1:n, k)
        subseq = perm[inds]
        pattern = Tuple(standardize(subseq))
        idx = perm_index[pattern]
        freq[idx] += 1
    end

    return freq
end

function profile_bruteforce_batched(batched_perms, order)
    m = length(batched_perms)
    A = Vector{Any}(undef, m)    
    @threads for i in 1:m
        p = batched_perms[i]
        A[i] = profile_bruteforce(p, order)
    end
    return A
end

