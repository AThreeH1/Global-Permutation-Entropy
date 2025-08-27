using JLD2
using Combinatorics
using Base.Threads
using BenchmarkTools
using Random
include("counting_functions.jl")

# Load JLD2 data instead of pickle
@load "./profile_dicts/profile_2.jld2" julia_profiles
profile_2_loaded = julia_profiles

all_keys_upto_2 = Set()
for lc in profile_2_loaded
    union!(all_keys_upto_2, keys(lc))
end

function profile_2(perm::Vector{Int})
    occurrences_fast = Dict{Tuple{Any,Any}, Int}()
    for key in all_keys_upto_2
        _, key1 = key
        occurrences_fast[key] = sum(vertex(perm,key1))
    end
    occurrences  = Dict{Tuple{Int,Vararg{Int}}, Rational{Int}}()
    python_perms = collect(permutations(0:1))
    julia_perms = [tuple((p .+ 1)...) for p in python_perms] 
    for (t,sigma) in enumerate(julia_perms)
        vector = profile_2_loaded[t]
        occurrences[Tuple(sigma)] = 0
        for key in keys(vector)
            val = vector[key]
            occurrences[sigma] += occurrences_fast[key]*val
        end
    end

    return [ Int(occurrences[σ]) for σ in julia_perms ]

end

function profile_2_batched(perms::Vector{Vector{Int}})
    results = Vector{Vector{Int}}(undef, length(perms))
    
    Threads.@threads for i in 1:length(perms)
        results[i] = profile_2(perms[i])
    end
    
    return results
end


