using Combinatorics, Random
using BenchmarkTools
using Base.Threads

include("profile2.jl")
include("profile3.jl")
include("profile4.jl")
include("profile5.jl")
include("profile6.jl")
include("profiles_bruteforce.jl")

function sliding_windows(seq::Vector{Int}, window_size::Int)
    n = length(seq)
    if window_size > n
        error("Window size cannot be larger than sequence length")
    end
    return [begin
        window = seq[i:i+window_size-1]
        indices = collect(1:window_size)
        sort!(indices, by = j -> window[j])
        perm = Vector{Int}(undef, window_size)
        for rank in 1:window_size
            perm[indices[rank]] = rank
        end
        perm
    end for i in 1:(n - window_size + 1)]
end

function test_profiles()
    for n in 6:20
    # Random.seed!(42)  # For reproducibility
        permutation = randperm(n) 

        println("Testing profiles for permutation: ", permutation)
        profile_two = profile_2(permutation)
        @assert profile_two == profile_bruteforce(permutation, 2)

        profile_three = profile_3(permutation)
        @assert profile_three == profile_bruteforce(permutation, 3)

        profile_four = profile_4(permutation)
        @assert profile_four == profile_bruteforce(permutation, 4)

        profile_five = profile_5(permutation)
        @assert profile_five == profile_bruteforce(permutation, 5)

        profile_six = profile_6(permutation)
        @assert profile_six == profile_bruteforce(permutation, 6)

    end

    # println("Timing profiles for permutation: ", permutation)
    permutation = randperm(20)

    println("Benchmarking profiles...")
    print("Time for profile_2: ")
    @btime profile_2($permutation)
    print("Time for profile_3: ")
    @btime profile_3($permutation)
    print("Time for profile_4: ")
    @btime profile_4($permutation)
    print("Time for profile_5: ")
    @btime profile_5($permutation)
    print("Time for profile_6: ")
    @btime profile_6($permutation)

    print("Time for profile_5_parallel: ")
    @btime profile_5_parallel($permutation)
    print("Time for profile_6_parallel: ")
    @btime profile_6_parallel($permutation)

    perm = randperm(100)
    batch_perms = sliding_windows(perm, 30)
    println("size of batch_perms: ", length(batch_perms), " of length ", length(batch_perms[1]))
    print("Time for batched profile 2:")
    @btime profile_2_batched($batch_perms)
    print("Time for batched profile 3:")
    @btime profile_3_batched($batch_perms)
    print("Time for batched profile 4:")
    @btime profile_4_batched($batch_perms)

end

if abspath(PROGRAM_FILE) == @__FILE__
    test_profiles()
end

test_profiles()