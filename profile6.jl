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

# Load JLD2 data instead of pickle
@load "./profile_dicts/profile_6.jld2" julia_profiles
profile_6_quadratic_loaded = julia_profiles

all_keys_upto_6_2 = Set()
for lc in profile_6_quadratic_loaded
    union!(all_keys_upto_6_2, keys(lc))
end

function profile_6(perm::Vector{Int};
                    only_6::Bool = true)
  
    occurrences_fast = Dict{Tuple{Any,Any}, Int}()
    for key in all_keys_upto_6_2
        key0, key1 = key
        if key0=="tree algo" 
            data = parse_python_literal(key1)
            mode, d1, d2, d3, d4, d5 = data
            perm0 = 
                mode == "e"  ? perm  :
                mode == "R"  ? rev1(perm) :
                mode == "C"  ? rev2(perm) :
                mode == "RC" ? rev1(rev2(perm)) :
                error("Unknown mode $mode")
            occurrences_fast[key] = count_gen(perm0, d1, d2, d3, d4, d5)

        elseif key0=="tree algo no middle"
            data = parse_python_literal(key1)
            mode, d1, d2, d3, d4, d5 = data
            perm0 = 
                mode == "e"  ? perm  :
                mode == "R"  ? rev1(perm) :
                mode == "C"  ? rev2(perm) :
                mode == "RC" ? rev1(rev2(perm)) :
                error("Unknown mode $mode")

            occurrences_fast[key] = count_gen(perm0, d1, d2, d3, d4, d5, no_middle=true)

        elseif key1=="count_NE_patterntree"
            if key0 == "e"
                occurrences_fast[key] = count_NE_patterntree(perm)
            elseif key0 == "R"
                occurrences_fast[key] = count_NE_patterntree(rev1(perm))
            elseif key0 == "C"
                occurrences_fast[key] = count_NE_patterntree(rev2(perm))
            end
    
        elseif key1=="count_SW_patterntree"
            if key0 == "e"
                occurrences_fast[key] = count_SW_patterntree(perm)
            elseif key0 == "R"
                occurrences_fast[key] = count_SW_patterntree(rev1(perm))
            elseif key0 == "C"
                occurrences_fast[key] = count_SW_patterntree(rev2(perm))
            end

        elseif key1=="marked_count_43215"
            if key0 == "e"
                occurrences_fast[key] = marked_count_43215(perm)
            elseif key0 == "R"
                occurrences_fast[key] = marked_count_43215(rev1(perm))
            elseif key0 == "C"
                occurrences_fast[key] = marked_count_43215(rev2(perm))
            elseif key0 == "RC"
                occurrences_fast[key] = marked_count_43215(rev1(rev2(perm)))
            end

        elseif key0=="pattern tree six"
            
            pattern_tree1 = [collect(y) for y in key1]
            # println("pattern_tree1 = ", pattern_tree1)
            occurrences_fast[key] = pt_count_level_6_quadratic(perm,pattern_tree1)

        else
           
            occurrences_fast[key] = sum(vertex(perm,key1))

        end
    end
   
    occurrences  = Dict{Tuple{Int,Vararg{Int}}, Rational{Int}}()

    if only_6
        python_perms = collect(permutations(0:5))
        julia_perms = [tuple((p .+ 1)...) for p in python_perms] 
        for (t ,sigma) in enumerate(julia_perms)
            vector = profile_6_quadratic_loaded[t+153]
            sigma_key = Tuple(sigma)
            occurrences[sigma_key] = 0
            for key in keys(vector)
                val = vector[key]
                occurrences[sigma] += occurrences_fast[key]*val
            end
        end

    else
        allperms_0based = Iterators.flatten(permutations(0:(k-1)) for k in 1:6)
        allperms_1based = [tuple((collect(p) .+ 1)...) for p in allperms_0based]
        for (t,sigma) in enumerate(allperms_1based)
            vector = profile_6_quadratic_loaded[t-1]
            sigma_key = Tuple(sigma)
            occurrences[sigma_key] = 0
            for key in keys(vector)
                val = vector[key]
                occurrences[sigma] += occurrences_fast[key]*val
            end
        end
    end
    
    return [ Int(occurrences[sigma]) for sigma in julia_perms]

end

function preprocess_keys_by_type_6()
    tree_algo_payloads = Any[]
    tree_algo_no_middle_payloads = Any[]
    count_NE_patterntree_payloads = Any[]
    count_SW_patterntree_payloads = Any[]
    marked_count_43215_payloads = Any[]
    pattern_tree_six_payloads = Any[]
    corner_tree_payloads = Any[]
    
    for key in all_keys_upto_6_2
        tag, payload = key
        if tag == "tree algo"
            push!(tree_algo_payloads, payload)
        elseif tag == "tree algo no middle"
            push!(tree_algo_no_middle_payloads, payload)
        elseif payload == "count_NE_patterntree"
            push!(count_NE_patterntree_payloads, tag)
        elseif payload == "count_SW_patterntree"
            push!(count_SW_patterntree_payloads, tag)
        elseif payload == "marked_count_43215"
            push!(marked_count_43215_payloads, tag)
        elseif tag == "pattern tree six"
            push!(pattern_tree_six_payloads, payload)
        else
            push!(corner_tree_payloads, payload)
        end
    end
    
    return tree_algo_payloads, tree_algo_no_middle_payloads, count_NE_patterntree_payloads, count_SW_patterntree_payloads, marked_count_43215_payloads, pattern_tree_six_payloads, corner_tree_payloads
end

TREE_ALGO_6_PAYLOADS, TREE_ALGO_NO_MIDDLE_6_PAYLOADS, COUNT_NE_6_PAYLOADS, COUNT_SW_6_PAYLOADS, MARKED_43215_6_PAYLOADS, PATTERN_TREE_SIX_6_PAYLOADS, CORNER_TREE_6_PAYLOADS = preprocess_keys_by_type_6()

function profile_6_parallel(perm::Vector{Int};
                            only_6::Bool = true, to_dict::Bool = false)
    
    # Give each thread its own copy of perm to avoid cache thrashing
    n_threads = Threads.nthreads()
    thread_perms = [copy(perm) for _ in 1:n_threads]
    
    # Thread-local dictionaries for each key type
    tree_algo_dicts = [Dict{Tuple{Any,Any}, Int}() for _ in 1:n_threads]
    tree_algo_no_middle_dicts = [Dict{Tuple{Any,Any}, Int}() for _ in 1:n_threads]
    count_NE_dicts = [Dict{Tuple{Any,Any}, Int}() for _ in 1:n_threads]
    count_SW_dicts = [Dict{Tuple{Any,Any}, Int}() for _ in 1:n_threads]
    marked_43215_dicts = [Dict{Tuple{Any,Any}, Int}() for _ in 1:n_threads]
    pattern_tree_six_dicts = [Dict{Tuple{Any,Any}, Int}() for _ in 1:n_threads]
    corner_tree_dicts = [Dict{Tuple{Any,Any}, Int}() for _ in 1:n_threads]
    
    # Process tree algo keys
    Threads.@threads for i in 1:length(TREE_ALGO_6_PAYLOADS)
        tid = Threads.threadid()
        payload = TREE_ALGO_6_PAYLOADS[i]
        
        data = parse_python_literal(payload)
        mode, d1, d2, d3, d4, d5 = data
        
        perm0 = if mode == "e"
            thread_perms[tid]
        elseif mode == "R"
            rev1(thread_perms[tid])
        elseif mode == "C"
            rev2(thread_perms[tid])
        elseif mode == "RC"
            rev1(rev2(thread_perms[tid]))
        else
            error("Unknown mode $mode")
        end
        
        tree_algo_dicts[tid][("tree algo", string(payload))] = count_gen(perm0, d1, d2, d3, d4, d5)
    end
    
    # Process tree algo no middle keys
    Threads.@threads for i in 1:length(TREE_ALGO_NO_MIDDLE_6_PAYLOADS)
        tid = Threads.threadid()
        payload = TREE_ALGO_NO_MIDDLE_6_PAYLOADS[i]
        
        data = parse_python_literal(payload)
        mode, d1, d2, d3, d4, d5 = data
        
        perm0 = if mode == "e"
            thread_perms[tid]
        elseif mode == "R"
            rev1(thread_perms[tid])
        elseif mode == "C"
            rev2(thread_perms[tid])
        elseif mode == "RC"
            rev1(rev2(thread_perms[tid]))
        else
            error("Unknown mode $mode")
        end
        
        tree_algo_no_middle_dicts[tid][("tree algo no middle", string(payload))] = count_gen(perm0, d1, d2, d3, d4, d5, no_middle=true)
    end
    
    # Process pattern tree five keys
    # Threads.@threads for i in 1:length(PATTERN_TREE_FIVE_6_PAYLOADS)
    #     tid = Threads.threadid()
    #     payload = PATTERN_TREE_FIVE_6_PAYLOADS[i]
        
    #     pattern_tree = [collect(sub) for sub in payload]
    #     pattern_tree_five_dicts[tid][("pattern tree five", payload)] = pt_count_level_5_quadratic(thread_perms[tid], pattern_tree)
    # end

    # Process count_NE_patterntree keys
    Threads.@threads for i in 1:length(COUNT_NE_6_PAYLOADS)
        tid = Threads.threadid()
        tag = COUNT_NE_6_PAYLOADS[i]

        perm0 = if tag == "e"
            thread_perms[tid]
        elseif tag == "R"
            reverse(thread_perms[tid])
        elseif tag == "C"
            [(length(thread_perms[tid])) - (x-1) for x in thread_perms[tid]]
        else
            error("Unknown tag $tag for count_NE_patterntree")
        end
        
        count_NE_dicts[tid][(tag, "count_NE_patterntree")] = count_NE_patterntree(perm0)
    end

    # Process count_SW_patterntree keys
    Threads.@threads for i in 1:length(COUNT_SW_6_PAYLOADS)
        tid = Threads.threadid()
        tag = COUNT_SW_6_PAYLOADS[i]

        perm0 = if tag == "e"
            thread_perms[tid]
        elseif tag == "R"
            reverse(thread_perms[tid])
        elseif tag == "C"
            [(length(thread_perms[tid])) - (x-1) for x in thread_perms[tid]]
        else
            error("Unknown tag $tag for count_SW_patterntree")
        end
        
        count_SW_dicts[tid][(tag, "count_SW_patterntree")] = count_SW_patterntree(perm0)
    end

    # Process marked_count_43215 keys
    Threads.@threads for i in 1:length(MARKED_43215_6_PAYLOADS)
        tid = Threads.threadid()
        tag = MARKED_43215_6_PAYLOADS[i]

        perm0 = if tag == "e"
            thread_perms[tid]
        elseif tag == "R"
            reverse(thread_perms[tid])
        elseif tag == "C"
            [(length(thread_perms[tid])) - (x-1) for x in thread_perms[tid]]
        elseif tag == "RC"
            reverse([(length(thread_perms[tid])) - (x-1) for x in thread_perms[tid]])
        else
            error("Unknown tag $tag for marked_count_43215")
        end
        
        marked_43215_dicts[tid][(tag, "marked_count_43215")] = marked_count_43215(perm0)
    end
    
    # Process pattern tree six keys
    Threads.@threads for i in 1:length(PATTERN_TREE_SIX_6_PAYLOADS)
        tid = Threads.threadid()
        payload = PATTERN_TREE_SIX_6_PAYLOADS[i]
        
        pattern_tree = [collect(sub) for sub in payload]
        pattern_tree_six_dicts[tid][("pattern tree six", payload)] = pt_count_level_6_quadratic(thread_perms[tid], pattern_tree)
    end
    
    # Process corner tree keys
    Threads.@threads for i in 1:length(CORNER_TREE_6_PAYLOADS)
        tid = Threads.threadid()
        payload = CORNER_TREE_6_PAYLOADS[i]
        
        corner_tree_dicts[tid][("corner tree", payload)] = sum(vertex(thread_perms[tid], payload))
    end
    
    # Merge all thread-local dictionaries
    occurrences_fast = Dict{Tuple{Any,Any}, Int}()
    
    for d in tree_algo_dicts
        for (k, v) in d
            occurrences_fast[k] = v
        end
    end
    
    for d in tree_algo_no_middle_dicts
        for (k, v) in d
            occurrences_fast[k] = v
        end
    end
    
    for d in count_NE_dicts
        for (k, v) in d
            occurrences_fast[k] = v
        end
    end

    for d in count_SW_dicts
        for (k, v) in d
            occurrences_fast[k] = v
        end
    end

    for d in marked_43215_dicts
        for (k, v) in d
            occurrences_fast[k] = v
        end
    end
    
    for d in pattern_tree_six_dicts
        for (k, v) in d
            occurrences_fast[k] = v
        end
    end
    
    for d in corner_tree_dicts
        for (k, v) in d
            occurrences_fast[k] = v
        end
    end

    # Rest of the computation (unchanged from profile_6)
    occurrences = Dict{Tuple{Int,Vararg{Int}}, Rational{Int}}()

    if only_6
        python_perms = collect(permutations(0:5))
        julia_perms = [tuple((p .+ 1)...) for p in python_perms] 
        for (t, σ) in enumerate(julia_perms)       
            vec = profile_6_quadratic_loaded[t + 153]
            σkey = Tuple(σ)
            occurrences[σkey] = 0//1
            for k in keys(vec)
                val = vec[k]
                occurrences[σkey] += occurrences_fast[k] * val
            end
        end
    else
        allperms_0based = Iterators.flatten(permutations(0:(k-1)) for k in 1:6)
        allperms_1based = [tuple((collect(p) .+ 1)...) for p in allperms_0based]
        for (t, σ) in enumerate(allperms_1based)
            vec = profile_6_quadratic_loaded[t-1]
            σkey = Tuple(σ)
            occurrences[σkey] = 0//1
            for k in keys(vec)
                val = vec[k]
                occurrences[σkey] += occurrences_fast[k] * val
            end
        end
    end

    return [round(Int, occurrences[σ]) for σ in julia_perms]
end

# A = [ 9,  0, 10,  8,  7,  1,  2, 13, 12, 14,  5, 11,  4,  3,  6] .+ 1
# B = profile_6(A)
# C = profile_6_parallel(A)
# println(B)
# println(C)
# println(sum(B))
# println(sum(C))
