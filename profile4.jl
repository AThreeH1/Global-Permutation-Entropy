using JLD2
using Combinatorics
using Base.Threads
using BenchmarkTools
using Random
include("counting_functions.jl")

function parse_python_literal(s::AbstractString)
    # Remove whitespace
    s = strip(s)
    
    # Handle tuples: (item1, item2, ...)
    if startswith(s, "(") && endswith(s, ")")
        inner = s[2:end-1]
        if isempty(strip(inner))
            return ()
        end
        
        # Split by commas, but be careful of nested structures
        items = split_python_items(inner)
        parsed_items = [parse_python_literal(strip(item)) for item in items]
        return tuple(parsed_items...)
    end
    
    # Handle lists: [item1, item2, ...]
    if startswith(s, "[") && endswith(s, "]")
        inner = s[2:end-1]
        if isempty(strip(inner))
            return []
        end
        
        items = split_python_items(inner)
        return [parse_python_literal(strip(item)) for item in items]
    end
    
    # Handle strings: 'string' or "string"
    if (startswith(s, "'") && endswith(s, "'")) || (startswith(s, "\"") && endswith(s, "\""))
        # Remove quotes and handle basic escape sequences
        inner = s[2:end-1]
        # Basic escape sequence handling (add more as needed)
        inner = replace(inner, "\\n" => "\n")
        inner = replace(inner, "\\t" => "\t")
        inner = replace(inner, "\\r" => "\r")
        inner = replace(inner, "\\'" => "'")
        inner = replace(inner, "\\\"" => "\"")
        inner = replace(inner, "\\\\" => "\\")
        return inner
    end
    
    # Handle numbers
    if occursin(r"^[+-]?\d+$", s)
        return parse(Int, s)
    end
    
    if occursin(r"^[+-]?\d*\.\d+([eE][+-]?\d+)?$", s)
        return parse(Float64, s)
    end
    
    # Handle booleans
    if s == "True"
        return true
    elseif s == "False"
        return false
    elseif s == "None"
        return nothing
    end
    
    # If nothing matches, return as string
    return s
end

# Helper function to split items while respecting nested brackets/parentheses
function split_python_items(s::AbstractString)
    items = String[]
    current_item = ""
    depth = 0
    in_string = false
    string_char = '\0'
    
    i = 1
    while i <= length(s)
        char = s[i]
        
        if !in_string && (char == '\'' || char == '"')
            in_string = true
            string_char = char
            current_item *= char
        elseif in_string && char == string_char
            # Check if it's escaped
            if i > 1 && s[i-1] == '\\'
                current_item *= char
            else
                in_string = false
                current_item *= char
            end
        elseif !in_string && (char == '(' || char == '[')
            depth += 1
            current_item *= char
        elseif !in_string && (char == ')' || char == ']')
            depth -= 1
            current_item *= char
        elseif !in_string && char == ',' && depth == 0
            push!(items, current_item)
            current_item = ""
        else
            current_item *= char
        end
        i += 1
    end
    
    if !isempty(current_item)
        push!(items, current_item)
    end
    
    return items
end

@load "./profile_dicts/profile_4.jld2" julia_profiles
profile_4_loaded = julia_profiles

all_keys_upto_4 = Set()
for lc in profile_4_loaded
    union!(all_keys_upto_4, keys(lc))
end

function profile_4(perm::Vector{Int})
    occurrences_fast = Dict{Tuple{Any,Any}, Int}()
    for key in all_keys_upto_4
        key0, key1 = key
        if key0=="corner tree"
            occurrences_fast[key] = sum(vertex(perm, key1))
        else
            data = parse_python_literal(key1)
            occurrences_fast[key] = count_gen(perm,data[2],data[3],data[4],data[5],data[6])
        end
    end
    occurrences  = Dict{Tuple{Int,Vararg{Int}}, Rational{Int}}()
    python_perms = collect(permutations(0:3))
    julia_perms = [tuple((p .+ 1)...) for p in python_perms] 
    for (t,sigma) in enumerate(julia_perms)
        vector = profile_4_loaded[t+8]
        occurrences[Tuple(sigma)] = 0//1
        for key in keys(vector)
            val = vector[key]
            occurrences[sigma] += occurrences_fast[key]*val
        end
    end

    return [ Int(occurrences[σ]) for σ in julia_perms ]

end

function preprocess_keys_by_type()
    tree_algo_payloads = Any[]
    corner_tree_payloads = Any[]
    
    for key in all_keys_upto_4
        tag, payload = key
        if tag == "tree algo"
            parsed_payload = parse_python_literal(payload)
            push!(tree_algo_payloads, (parsed_payload, string(payload)))
        else
            push!(corner_tree_payloads, payload)
        end
    end
    
    return tree_algo_payloads, corner_tree_payloads
end

# Call this once when your module loads
const TREE_ALGO_PAYLOADS, CORNER_TREE_PAYLOADS = preprocess_keys_by_type()

function profile_4_parallel(perm::Vector{Int})

    n_threads = Threads.nthreads()
    thread_perms = [copy(perm) for _ in 1:n_threads]

    tree_algo_dicts = [Dict{Tuple{Any,Any}, Int}() for _ in 1:n_threads]
    corner_tree_dicts = [Dict{Tuple{Any,Any}, Int}() for _ in 1:n_threads]
    
    # println("Total threads: ", Threads.nthreads())
    Threads.@threads for i in 1:length(CORNER_TREE_PAYLOADS)
        # println("Current thread: ", Threads.threadid())
        tid = Threads.threadid()
        payload = CORNER_TREE_PAYLOADS[i]
        corner_tree_dicts[tid][("corner tree", payload)] = sum(vertex(thread_perms[tid], payload))
    end

    Threads.@threads for i in 1:length(TREE_ALGO_PAYLOADS)
        # println("Current thread: ", Threads.threadid())
        tid = Threads.threadid()
        data, payload = TREE_ALGO_PAYLOADS[i]
        # data = parse_python_literal(payload)
        tree_algo_dicts[tid][("tree algo", payload)] = count_gen(thread_perms[tid], data[2], data[3], data[4], data[5], data[6])
    end
        
    # Merge thread-local dictionaries
    occurrences_fast = Dict{Tuple{Any,Any}, Int}()
    for d in corner_tree_dicts
        for (k, v) in d
            occurrences_fast[k] = v
        end
    end

    for d in tree_algo_dicts
        for (k, v) in d
            occurrences_fast[k] = v
        end
    end

    # Rest of the computation (unchanged from original)
    occurrences = Dict{Tuple{Int,Vararg{Int}}, Rational{Int}}()
    python_perms = collect(permutations(0:3))
    julia_perms = [tuple((p .+ 1)...) for p in python_perms] 
    
    for (t, sigma) in enumerate(julia_perms)
        vector = profile_4_loaded[t+8]
        occurrences[sigma] = 0//1
        for key in keys(vector)
            val = vector[key]
            occurrences[sigma] += occurrences_fast[key] * val
        end
    end

    return [Int(occurrences[σ]) for σ in julia_perms]
end

function profile_4_batched(perms::Vector{Vector{Int}})
    results = Vector{Vector{Int}}(undef, length(perms))
    
    Threads.@threads for i in 1:length(perms)
        results[i] = profile_4(perms[i])
    end
    
    return results
end

function main()

    seq_lengths = [50, 100, 200, 500, 1000]

    println("Threads available: ", nthreads())

    println("\n=== Internal Parallel vs Sequential Timing ===")
    for len in seq_lengths
        println("\n-- Sequence Length: $len --")
        Random.seed!(1234)
        perm = randperm(len)

        print("Sequential: ")
        @btime profile_4($perm)

        print("Internal Parallel: ")
        @btime profile_4_parallel($perm)

    end

    batch_sizes = [4, 8, 12, 16, 32]

    for batch in batch_sizes
        println("\n-- Batch Count: $batch --")
        Random.seed!(1234)
        perms = [randperm(200) for _ in 1:batch]

        print("Time Taken for batched profile 2, Seq Len 500: ")
        @btime profile_4_batched($perms)
    end
end

# main()
