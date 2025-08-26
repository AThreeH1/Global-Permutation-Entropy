using Combinatorics, Random
using BenchmarkTools
using Base.Threads

include("profile2.jl")
include("profile3.jl")
include("profile4.jl")
include("profile5.jl")
include("profile6.jl")
include("profiles_bruteforce.jl")

function count_brute_force(pattern, permutation)
    count = 0
    for sigma in combinations(permutation, length(pattern))
        if collect(pattern) == sortperm(sigma)
            count += 1
        end
    end
    return count
end

function profile_level_two_brute_force(permutation)
    occurrences = Dict{Vector{Int}, Int}()
    for sigma in permutations(1:2)
        occurrences[collect(sigma)] = count_brute_force(sigma, permutation)
    end
    result = [occurrences[s] for s in permutations(1:2)]
    return [Int(x) for x in result]
end

function profile_level_three_brute_force(permutation)
    occurrences = Dict{Tuple{Int, Int, Int}, Int}()
    python_perms = collect(permutations(0:2))
    julia_perms = [tuple((p .+ 1)...) for p in python_perms]
    for sigma in julia_perms
        occurrences[sigma] = count_brute_force(sigma, permutation)
    end
    result = [occurrences[s] for s in julia_perms]
    return [Int(x) for x in result]
end

function profile_level_four_brute_force(permutation)
    occurrences = Dict{Tuple{Int, Int, Int, Int}, Int}()
    python_perms = collect(permutations(0:3))
    julia_perms = [tuple((p .+ 1)...) for p in python_perms] 
    for sigma in julia_perms
        occurrences[sigma] = count_brute_force(sigma, permutation)
    end
    result = [occurrences[s] for s in julia_perms]
    return [Int(x) for x in result]
end

function profile_level_five_brute_force(permutation)
    occurrences = Dict{Tuple{Int, Int, Int, Int, Int}, Int}()
    python_perms = collect(permutations(0:4))
    julia_perms = [tuple((p .+ 1)...) for p in python_perms] 
    for sigma in julia_perms
        occurrences[sigma] = count_brute_force(sigma, permutation)
    end
    result = [occurrences[s] for s in julia_perms]
    return [Int(x) for x in result]
end

function profile_level_six_brute_force(permutation)
    occurrences = Dict{Tuple{Int, Int, Int, Int, Int, Int}, Int}()
    python_perms = collect(permutations(0:5))
    julia_perms = [tuple((p .+ 1)...) for p in python_perms] 
    for sigma in julia_perms
        occurrences[sigma] = count_brute_force(sigma, permutation)
    end
    result = [occurrences[s] for s in julia_perms]
    return [Int(x) for x in result]
end

function batch_profile_level_six(perms::Vector{Vector{Int}})
    B = length(perms)
    results = Vector{Vector{<:Any}}(undef, B)  
    @threads for i in 1:B
        p = perms[i]
        res = profile_level_six_with_dp_quadr(p)
        results[i] = res
    end
    return results
end

function gen_permutations(n::Int, batch::Int)
    [randperm(n) for _ in 1:batch]
end

function test_profiles()
    for n in 5:50
    # Random.seed!(42)  # For reproducibility
        permutation = randperm(n) 

        println("Testing profiles for permutation: ", permutation)
        # profile_two = sort(profile_2(permutation))
        # @assert profile_two == sort(profile_bruteforce(permutation, 2))

        # profile_three = sort(profile_3(permutation))
        # @assert profile_three == sort(profile_level_three_brute_force(permutation))

        # profile_four = sort(profile_4(permutation))
        # @assert profile_four == sort(profile_level_four_brute_force(permutation))

        profile_five = profile_5(permutation)
        @assert profile_five == profile_bruteforce(permutation, 5)

        # profile_six = sort(profile_6(permutation))
        # @assert profile_six == sort(profile_level_six_brute_force(permutation))

    end

    # println("Timing profiles for permutation: ", permutation)
    # @btime profile_2($permutation)
    # @btime profile_3($permutation)
    # @btime profile_4($permutation)
    # @btime profile_5($permutation)
    # @btime profile_6($permutation)
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_profiles()
end

test_profiles()