include("profile2.jl")
include("profile3.jl")
include("profile4.jl")
include("profile5.jl")
include("profile6.jl")
include("profiles_bruteforce.jl")

function global_permutation_entropy(sequence::Vector{<:Real}; order::Int, window_size::Int64, 
                                   stride::Int64 = 1, normalized::Bool = true, parallel::Bool = false, bruteforce::Bool=false)
    
    if bruteforce==true

        n = length(sequence)
        
        # Check if window_size is valid
        if window_size > n
            error("Window size cannot be larger than sequence length")
        end
        
        # Calculate number of windows
        num_windows = div(n - window_size, stride) + 1
        entropies = Vector{Float64}(undef, num_windows)
        
        # Create batches of windows
        windows = Vector{Vector{Int}}(undef, num_windows)
        
        for i in 1:num_windows
            start_idx = (i - 1) * stride + 1
            end_idx = start_idx + window_size - 1
            window = sequence[start_idx:end_idx]
            
            # Convert window to permutation representation
            indices = collect(1:window_size)
            sort!(indices, by = j -> window[j])
            perm = Vector{Int}(undef, window_size)
            for rank in 1:window_size
                perm[indices[rank]] = rank
            end
            
            windows[i] = perm
        end
        
        occurrences_batch = profile_bruteforce_batched(windows, order)

        # Calculate entropy for each window
        for i in 1:num_windows
            occurrences = occurrences_batch[i]
            total = sum(occurrences)
            probs = occurrences ./ total
            entropy = 0.0
            
            for p in probs
                if p > 0
                    entropy -= p * log2(p)
                end
            end
            
            if normalized
                log2_num_permutations = log2(factorial(order))
                entropies[i] = entropy / log2_num_permutations
            else
                entropies[i] = entropy
            end
        end

    else

        n = length(sequence)
        
        # Check if window_size is valid
        if window_size > n
            error("Window size cannot be larger than sequence length")
        end
        
        # Calculate number of windows
        num_windows = div(n - window_size, stride) + 1
        entropies = Vector{Float64}(undef, num_windows)
        # If order is 2, 3, or 4, we can use batched profiles
        if order == 2 || order == 3 || order == 4
            # Create batches of windows
            windows = Vector{Vector{Int}}(undef, num_windows)
            
            for i in 1:num_windows
                start_idx = (i - 1) * stride + 1
                end_idx = start_idx + window_size - 1
                window = sequence[start_idx:end_idx]
                
                # Convert window to permutation representation
                indices = collect(1:window_size)
                sort!(indices, by = j -> window[j])
                perm = Vector{Int}(undef, window_size)
                for rank in 1:window_size
                    perm[indices[rank]] = rank
                end
                
                windows[i] = perm
            end
            
            if order == 2
                occurrences_batch = profile_2_batched(windows)
            elseif order == 3
                occurrences_batch = profile_3_batched(windows)
            elseif order == 4
                occurrences_batch = profile_4_batched(windows)
            else
                error("Order must be between 2 and 4")
            end
            
            # Calculate entropy for each window
            for i in 1:num_windows
                occurrences = occurrences_batch[i]
                total = sum(occurrences)
                probs = occurrences ./ total
                entropy = 0.0
                
                for p in probs
                    if p > 0
                        entropy -= p * log2(p)
                    end
                end
                
                if normalized
                    log2_num_permutations = log2(factorial(order))
                    entropies[i] = entropy / log2_num_permutations
                else
                    entropies[i] = entropy
                end
            end
            

        else
        # If order is 5 and 6, we can use internally parallelizes profiles
            
            for i in 1:num_windows
                start_idx = (i - 1) * stride + 1
                end_idx = start_idx + window_size - 1
                window = sequence[start_idx:end_idx]
                
                # Convert window to permutation representation
                indices = collect(1:window_size)
                sort!(indices, by = j -> window[j])
                perm = Vector{Int}(undef, window_size)
                for rank in 1:window_size
                    perm[indices[rank]] = rank
                end
                
                # Get occurrences for this window
                if parallel == false
                    if order == 5
                        occurrences = profile_5(perm)
                    elseif order == 6
                        occurrences = profile_6(perm)
                    else
                        error("Order must be between 2 and 6")
                    end
                else
                    if order == 5
                        occurrences = profile_5_parallel(perm)
                    elseif order == 6
                        occurrences = profile_6_parallel(perm)
                    else
                        error("Order must be between 2 and 6 for parallel processing")
                    end
                end
                
                # Calculate entropy for this window
                total = sum(occurrences)
                probs = occurrences ./ total
                entropy = 0.0
                
                for p in probs
                    if p > 0
                        entropy -= p * log2(p)
                    end
                end
                
                if normalized
                    log2_num_permutations = log2(factorial(order))
                    entropies[i] = entropy / log2_num_permutations
                else
                    entropies[i] = entropy
                end
            end
        end
    end

    return entropies
end

function main()
    seq = rand(1:100, 50)
    # println(global_permutation_entropy(seq, 2, 10))
    # println(global_permutation_entropy(seq, 3, 20))
    # println(global_permutation_entropy(seq, 4, 30))
    println(global_permutation_entropy(seq, 5, 40))
    println(global_permutation_entropy(seq, 5, 40, bruteforce=true))
    # println(global_permutation_entropy(seq, 6, 50))
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("Running main function...")
    main()
end

main()
