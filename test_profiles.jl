using Combinatorics, Random
using BenchmarkTools
using Base.Threads

include("profile2.jl")
include("profile3.jl")
include("profile4.jl")
include("profile5.jl")
include("profile6.jl")
include("profiles_bruteforce.jl")




function test_profiles()
    for n in 6:30
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