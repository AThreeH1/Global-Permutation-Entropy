using Random
using Statistics
using Plots
using BenchmarkTools
using Base.Threads
const NE, NW, SE, SW = "NE", "NW", "SE", "SW"
using Combinatorics
using Symbolics


mutable struct SumTree
    # Segment tree for range sum queries.
    # We only need suffix and prefix queries.
    length::Int
    arrays::Vector{Vector{Int}}
    
    function SumTree(n::Int)
        num_levels = ceil(Int, log2(max(1, n))) + 1
        arrays = [zeros(Int, 2^i) for i in 0:(num_levels-1)]
        reverse!(arrays)
        new(n, arrays)
    end
end


function add!(tree::SumTree, i::Int, val::Int)
    """Add val to position i in the tree"""
    if i < 1 || i > tree.length
        throw(BoundsError(tree, i))
    end
    
    i_zero = i - 1  # Convert to 0-based indexing
    for x in 1:length(tree.arrays)
        level_size = 2^(x-1)
        idx = div(i_zero, level_size) + 1
        if idx <= length(tree.arrays[x])
            tree.arrays[x][idx] += val
        end
    end
end


function sum_suffix(tree::SumTree, i::Int)
    """Sum of all elements at positions i+1 to n"""
    if i < 0 || i > tree.length
        throw(BoundsError(tree, i))
    end
    
    if i == tree.length
        return 0  # No elements after position n
    end
    
    i_zero = i - 1  # Convert to 0-based indexing
    total = 0
    for x in 1:(length(tree.arrays)-1)
        level_size = 2^(x-1)
        if div(i_zero, level_size) % 2 == 0 && div(i_zero, level_size) + 2 <= length(tree.arrays[x])
            total += tree.arrays[x][div(i_zero, level_size) + 2]
        end
    end
    return total
end


function sum_prefix(tree::SumTree, i::Int)
    """Sum of all elements at positions 1 to i-1"""
    if i < 1 || i > tree.length + 1
        throw(BoundsError(tree, i))
    end
    
    if i == 1
        return 0  # No elements before position 1
    end
    
    i_zero = i - 1  # Convert to 0-based indexing
    total = 0
    for x in 1:(length(tree.arrays)-1)
        level_size = 2^(x-1)
        if div(i_zero, level_size) % 2 == 1 && div(i_zero, level_size) >= 1
            total += tree.arrays[x][div(i_zero, level_size)]
        end
    end
    return total
end

"""
==============================================================================================================
counting occurences of corner trees: vertex and edge functions
==============================================================================================================
"""
function vertex(permutation::Vector{Int}, tree)
    n = length(permutation)
    A = ones(Int, n)
    
    for (edge_label, child) in tree
        child_result = edge(permutation, child, edge_label)
        A = A .* child_result
    end
    
    return A
end


function edge(permutation::Vector{Int}, tree, edge_label::String)
    n = length(permutation)
    A = vertex(permutation, tree)
    C = zeros(Int, n)
    B = SumTree(n)
    S, W = edge_label
    
    if W == 'W'
        for i in 1:n
            idx = permutation[i]
            if S == 'S'
                tmp = sum_prefix(B, idx)
                C[i] = tmp
            else
                C[i] = sum_suffix(B, idx)
            end
            add!(B, idx, A[i])
        end
    else
        for i in n:-1:1
            idx = permutation[i]
            if S == 'S'
                C[i] = sum_prefix(B, idx)
            else
                C[i] = sum_suffix(B, idx)
            end
            add!(B, idx, A[i])
        end
    end
    
    return C
end

"""
==============================================================================================================
counting occurences of pure west corner tree
==============================================================================================================
"""

function countW(permutation::Vector{Int}, tree;
    row::Union{Nothing,Int}=nothing,
    chunk_size::Union{Nothing,Int}=nothing,
    not_B::Bool=false)

    # — Shared state for the two inner functions —
    A_v_dict    = Dict{Any, Vector{Int}}()
    B_e_dict    = Dict{Any, SumTree}()
    touched_dict = Dict{Tuple{Any,Int}, Bool}()

    # Inner function: like your Python vertex_W
    function vertex_W(subtree, i::Int)
        # init if first time we see this tree
        if !haskey(A_v_dict, subtree)
            A_v_dict[subtree] = ones(Int, length(permutation))
        end
        # only compute once per (tree,i)
        if !get(touched_dict, (subtree, i), false)
            touched_dict[(subtree, i)] = true
            if isa(subtree, Tuple)
                for entry in subtree
                    if isa(entry, Tuple) && length(entry) == 2
                        label, child = entry
                        A_v_dict[subtree][i] *= edge_W(child, label, i)
                    end
                end
            end
        end
        return A_v_dict[subtree][i]
    end

    # Inner function: like your Python edge_W
    function edge_W(child, edge_label, i::Int)
        edge_id = (child, edge_label)

        # init SumTree on first visit
        if !haskey(B_e_dict, edge_id)
            B_e_dict[edge_id] = SumTree(length(permutation))
        end

        idx = permutation[i]
        if !get(touched_dict, (edge_id, i), false)
            touched_dict[(edge_id, i)] = true
            # assume SumTree has a mutating add! method
            add!(B_e_dict[edge_id], idx, vertex_W(child, i))
        end

        # check first character of the label
        return edge_label[1] == 'S' ?
            sum_prefix(B_e_dict[edge_id], idx) :
            sum_suffix(B_e_dict[edge_id], idx)
    end

    ret = Int[]  # result vector

    if not_B
        # Python: if not_B: …
        @assert row !== nothing && chunk_size !== nothing
        countCT = 0
        countCTback = 0
        for i in 1:length(permutation)
            y = permutation[i]
            # align to chunk boundaries (0-based logic → (i-1) % chunk)
            if (i-1) % chunk_size == 0
                countCTback = countCT
            end
            if y < row
                countCT += vertex_W(tree, i)
            end
            if row <= y && y < row + chunk_size
                push!(ret, countCT - countCTback)
            else
                push!(ret, 0)
            end
        end
        # println("ret1 = ", ret)

    elseif row !== nothing
        # Python: elif row is not None: …
        @assert chunk_size !== nothing
        countCTW = 0
        for i in 1:length(permutation)
            y = permutation[i]
            if y < row
                countCTW += vertex_W(tree, i)
            end
            if row <= y < row + chunk_size
                push!(ret, countCTW)
            else
                push!(ret, 0)
            end
        end

    else
        # Python: else: …
        for i in 1:length(permutation)
        push!(ret, vertex_W(tree, i))
        end
    end

    return ret
end

function count_Box(permutation::Vector{Int},
    left_trees,
    middle_tree,
    right_trees;
    no_middle::Bool=false)

    n     = length(permutation)
    chunk = floor(Int, n^(1/3))
    permT = invperm(permutation)    # assumes you have invperm defined

    count = 0

    # initialize the “boxes”
    right_boxes  = [NDProductTree(n, 2) for _ in right_trees]
    left_boxes   = [NDProductTree(n, 2) for _ in left_trees]

    middle_box   = no_middle ? nothing : NDProductTree(n, 2)

    # Fill left_boxes
    for (i, tree) in pairs(left_trees)
        outcome = countW(permutation, tree)     # pure mode
        for x in 1:n
            y = permutation[x]
            add_NDP(left_boxes[i], [x, y], outcome[x])
        end
    end

    # Fill right_boxes
    for (i, tree) in pairs(right_trees)
        outcome = countW(permutation, tree)     # pure mode
            for x in 1:n
                y = permutation[x]
                add_NDP(right_boxes[i], [x, y], outcome[x])
        end
    end

    # Fill middle_box if requested
    if !no_middle
        outcome = countW(permutation, middle_tree)
            for x in 1:n
                y = permutation[x]
                add_NDP(middle_box, [x, y], outcome[x])
        end
    end

    # Main counting loop
    for x4 in 1:n
        y4 = permutation[x4]
        # align to the chunk grid (1-based version)
        col = x4 - ((x4 - 1) % chunk)
        row = y4 - ((y4 - 1) % chunk)

        # x1 from col up to x4-1 (since Python slice col:x4 excluded x4)
        for x1 in col:(x4-1)
            y1 = permutation[x1]

            # y3 loop: indices row to y4-1
            for y3 in row:(y4-1)
                x3 = permT[y3]

                if x3 < x1 && y3 > y1
                    # accumulate products over left and right boxes
                    left_product = 1
                    right_product = 1
                    for box in left_boxes
                        left_product *= sum_box_NDP(box, [(0, x1), (0, y1)])
                    end 
                    for box in right_boxes
                        right_product *= sum_box_NDP(box, [(0, x3), (0, y3)])
                    end
                    if !no_middle
                        count += left_product * right_product *
                                sum_box_NDP(middle_box, [(x3+1, x1), (y1+1, y3)])
                    else
                        count += left_product * right_product
                    end
                end
            end
        end
    end

    return count
end




"""
==============================================================================================================
counting occurences of TreeAlgo_5/3
==============================================================================================================
"""
function count_gen(permutation::Vector{Int},
                   tree,
                   inverse_tree,
                   left_trees,
                   middle_tree,
                   right_trees;
                   no_middle::Bool=false)

    n     = length(permutation)
    chunk = floor(Int, n^(1/3))
    permT = invperm(permutation)

    total = 0

    # Sum over “pure” blocks for the main tree
    for row in 1:chunk:(n)
        total += sum(
            countW(permutation, tree;
                   row=row,
                   chunk_size=chunk,
                   not_B=false)
        )
    end

    # Sum over “not_B” blocks for the inverse_tree on permT
    for col in 1:chunk:(n)
        total += sum(
            countW(permT, inverse_tree;
                   row=col,
                   chunk_size=chunk,
                   not_B=true)
        )
    end

    # Finally, add the box‐based counts (with or without middle)
    total += count_Box(permutation,
                       left_trees,
                       middle_tree,
                       right_trees;
                       no_middle=no_middle)

    return total
end


function dyadic_range_1d(l::Int, r::Int)
    x = 1
    result = Tuple{Int, Int}[]
    while l < r
        if l % (2x) != 0
            push!(result, (l, l + x))
            l += x
        elseif r % (2x) != 0
            push!(result, (r - x, r))
            r -= x
        else
            x *= 2
        end
    end
    return result
end

function dyadic_ranges_nd(bounds::Vector{Tuple{Int,Int}})
    ranges1d = Vector{Vector{Tuple{Int,Int}}}(undef, length(bounds))
    @inbounds for i in 1:length(bounds)
        l, r = bounds[i]
        ranges1d[i] = collect(dyadic_range_1d(l, r))
    end
    return Iterators.product(ranges1d...)
end

"""
==============================================================================================================
For query interval sums efficiently using dyadics
==============================================================================================================
"""

mutable struct NDProductTree
    # Structure for efficient n-dimensional range/box sum queries.
    n::Int
    ndim::Int
    logn::Int
    table::Dict{Any,Int}

    function NDProductTree(n::Int, ndim::Int)
        logn  = ceil(Int, log2(max(1, n)))
        table = Dict{Any,Int}()  

        if ndim == 2
            sizehint!(table, n^2 * (logn + 1)^2)
        elseif ndim == 4
            sizehint!(table, n^2 * (logn + 1)^4)
        end
        new(n, ndim, logn, table)
    end
end


function add_NDP(tree::NDProductTree, coords::Vector{Int}, value::Int)
    @assert length(coords) == tree.ndim
    # time_start = time_ns()
    # build dyadic ranges along each axis
    ranges_per_dim = Vector{Vector{Tuple{Int,Int}}}(undef, tree.ndim)
    @inbounds for (i, x) in enumerate(coords)
        axis_ranges = Vector{Tuple{Int,Int}}(undef, tree.logn + 1)
        for level in 0:tree.logn
            size   = 1 << level 
            start  = x & ~(size - 1) 
            finish = start + size
            axis_ranges[level + 1] = (start, finish)
        end
        ranges_per_dim[i] = axis_ranges
    end

    if tree.ndim == 2
        @inbounds for r1 in ranges_per_dim[1]
            for r2 in ranges_per_dim[2]
                key = (r1, r2)
                tree.table[key] = get(tree.table, key, 0) + value
            end
        end
    elseif tree.ndim == 4
        @inbounds for r1 in ranges_per_dim[1]
            for r2 in ranges_per_dim[2]
                for r3 in ranges_per_dim[3]
                    for r4 in ranges_per_dim[4]
                        key = (r1, r2, r3, r4)
                        tree.table[key] = get(tree.table, key, 0) + value
                    end
                end
            end
        end
    else
        # Fallback for other dimensions
        for key in Iterators.product(ranges_per_dim...)
            kt = Tuple(key)
            tree.table[kt] = get(tree.table, kt, 0) + value
        end
    end
end



function sum_box_NDP(tree::NDProductTree, bounds)
    @assert length(bounds) == tree.ndim
    total = 0
    @inbounds for key in dyadic_ranges_nd(bounds)
        kt = Tuple(key)
        total += get(tree.table, kt, 0)
    end
    return total
end

function shape(c)
    sc = sort(c)
    return [findfirst(==(x), sc) for x in c]
end


function generate_matches(permutation, pattern)
    r = length(pattern)
    matches = []
    @inbounds for indices in combinations(1:length(permutation), r)
        values = permutation[indices]
        if shape(values) == pattern
            coords = collect(zip(indices .- 1, values))
            push!(matches, coords)
        end
    end
    return matches
end

@variables x0 x1 y0 y1 n

function evaluate_expr(symbolic_input, values::Dict{String, Int})
    expr_sym = if symbolic_input isa String
        ex = Meta.parse(symbolic_input)     
        eval(ex)                            
    else
        symbolic_input
    end

    vars = Symbolics.get_variables(expr_sym)  
    subs_map = Dict{Any,Any}()
    @inbounds for var in vars
        name = string(var)  
        if haskey(values, name)
            subs_map[var] = values[name]
        else
            error("No value provided for variable $name")
        end
    end

    substituted = Symbolics.substitute(expr_sym, subs_map)
    return Symbolics.value(substituted)  
end

"""
==============================================================================================================
For the counting of pattern trees with six vertices in quadratic time
==============================================================================================================
"""

function pt_count_level_6_quadratic(perm::Vector{Int}, pt_tree)
    n = length(perm) + 1
    pt_tree[1] = pt_tree[1] .+ 1
    pt_tree[2] = pt_tree[2] .+ 1
    pt_tree[3] = pt_tree[3] .+ 1
    # Tree for indices 2 and 3
    tree_two_three = NDProductTree(n, 4)
    pattern_two_three = pt_tree[2]
    matches_two_three = generate_matches(perm, pattern_two_three)
    @inbounds for match in matches_two_three
        (x2, y2), (x3, y3) = match
        add_NDP(tree_two_three, [x2, x3, y2, y3], 1)
    end

    # Tree for indices 4 and 5
    tree_four_five = NDProductTree(n, 4)
    pattern_four_five = pt_tree[3]
    @inbounds matches_four_five = generate_matches(perm, pattern_four_five)
    for match in matches_four_five
        (x4, y4), (x5, y5) = match
        add_NDP(tree_four_five, [x4, x5, y4, y5], 1)
    end

    # Count using patterns 0 and 1
    pattern_zero_one = pt_tree[1]
    matches_zero_one = generate_matches(perm, pattern_zero_one)
    count_DP_fast = 0

    @inbounds for match in matches_zero_one
        (x0, y0), (x1, y1) = match
        context = Dict(
            "n" => n,
            "x0" => x0,
            "x1" => x1,
            "y0" => y0,
            "y1" => y1
        )

        bounds_two_three = [
            (evaluate_expr(left, context), evaluate_expr(right, context))
            for (left, right) in pt_tree[4]
        ]

        bounds_four_five = [
            (evaluate_expr(left, context), evaluate_expr(right, context))
            for (left, right) in pt_tree[5]
        ]
        count_DP_fast += sum_box_NDP(tree_two_three, bounds_two_three) * sum_box_NDP(tree_four_five, bounds_four_five)
    end

    return count_DP_fast
end


function marked_count_3214(perm::Vector{Int})
    n = length(perm)
    chunk = Int(floor(n^(1/3)))
    permT = invperm(perm)
    T_out = NDProductTree(n, 2)
    # making perm and permT 0 index based
    perm = perm .-1
    permT = permT .-1

    #-------------------------------
    # 3,4 not in same row
    #-------------------------------

    for row in 0:chunk:n-1
        T1 = NDProductTree(n, 2)
        T2 = NDProductTree(n, 2)
        T3 = NDProductTree(n, 2)

        for i in 0:n-1
            if perm[i+1] < row
                add_NDP(T1, [i, perm[i+1]], 1)
                add_NDP(T2, [i, perm[i+1]], sum_box_NDP(T1, [(0, i), (perm[i+1], row)]))
                add_NDP(T3, [i, perm[i+1]], sum_box_NDP(T2, [(0, i), (perm[i+1], row)]))
            elseif row <= perm[i+1] < row + chunk
                add_NDP(T_out, [i, perm[i+1]], sum_box_NDP(T3, [(0, i), (0, row)]))
            end
        end
    end

    # -------------------------------------------------------------            
    # if 1,4 not in same col, but 3,4 in same row
    #--------------------------------------------------------------

    for col in (chunk):chunk:n-1
        TT1 = NDProductTree(n, 2)
        TT2 = NDProductTree(n, 2)
        TT3 = NDProductTree(n, 2)
        
        for (y, x) in enumerate(permT)
            if x < col
                add_NDP(TT1, [y-1, x], 1)
                add_NDP(TT2, [y-1, x], sum_box_NDP(TT1, [(0, y-1), (x, col)]))
                add_NDP(TT3, [y-1, x], sum_box_NDP(TT2, [(0, y-1), (x, col)]))
            elseif col <= x < col + chunk
                row = (y-1) - (y - 1) % chunk
                add_NDP(T_out, [x, y-1], sum_box_NDP(TT3, [(row, y-1), (0, col)]))
            end
        end
    end

    # ------------------------------------------------
    # Final section - both conditions
    # ------------------------------------------------

    PT = NDProductTree(n, 2)
    for (x2, y2) in enumerate(perm)
        add_NDP(PT, [x2-1, y2], 1)
    end
    
    for (x4, y4) in enumerate(perm)
        col = (x4-1) - (x4 - 1) % chunk
        row = y4 - y4 % chunk
        
        for x1 in col:(x4-2)
            y1 = perm[x1+1]
            for y3 in row:(y4-1)
                x3 = permT[y3+1]
                if x3 < x1 && y3 > y1
                    add_NDP(T_out, [x4-1, y4], sum_box_NDP(PT, [(x3+1, x1), (y1+1, y3)]))
                end
            end
        end
    end
    
    return T_out
end

function count_SW_patterntree(perm::Vector{Int})
    n = length(perm)
    T_sw = NDProductTree(n, 2)

    T_3214 = marked_count_3214(perm)

    perm = perm .-1

    for i in 0:n-1
        # Python: T_3214.sum_box([(0,i),(0,perm[i])])
        # Julia equivalent: [(1,i),(1,perm[i])]
        add_NDP(T_sw, [i, perm[i+1]], sum_box_NDP(T_3214, [(0, i), (0, perm[i+1])]))
    end
    
    # Python: T_sw.sum_box([(0,n),(0,n)])
    # Julia equivalent: [(1,n),(1,n)]
    return sum_box_NDP(T_sw, [(0, n), (0, n)])
end


"""
==============================================================================================================
Counts occurences of special pattern tree [1]---NE---[321(4)] in a permutation.
Tree Vector = 3*#3214 + #32145 + #42135 + #32415 + #41325 + #24315 + #14325 + 2*#42315 + 2*#34215 + 2*#43125 + 4*#43215
==============================================================================================================
"""


function count_NE_patterntree(perm::Vector{Int})
    n = length(perm)
    T_ne = NDProductTree(n, 2)

    T_3214 = marked_count_3214(perm)

    perm = perm .-1

    for i in 0:n-1
        add_NDP(T_ne, [i, perm[i+1]], sum_box_NDP(T_3214, [(i+1, n), (perm[i+1]+1, n)]))
    end
    
    return sum_box_NDP(T_ne, [(0, n), (0, n)])
end

# perm = [ 9,  0, 10,  8,  7,  1,  2, 13, 12, 14,  5, 11,  4,  3,  6] .+ 1
# println("this thing here = ", count_NE_patterntree(perm))


"""
=================================================================================================
A data structure that queries the number of ascending pairs within given bounds in O(q) time.
=================================================================================================
"""
mutable struct PairTree
    n::Int
    tree1::NDProductTree
    tree2::NDProductTree
    Vert::Dict{Int, Set{Tuple{Int,Int}}}
    Hori::Dict{Int, Set{Tuple{Int,Int}}}
    Vert_Trees::Dict{Int, NDProductTree}
    Hori_Trees::Dict{Int, NDProductTree}
    
    function PairTree(n::Int)
        new(n, 
            NDProductTree(n, 2), 
            NDProductTree(n, 2),
            Dict{Int, Set{Tuple{Int,Int}}}(),
            Dict{Int, Set{Tuple{Int,Int}}}(),
            Dict{Int, NDProductTree}(),
            Dict{Int, NDProductTree}())
    end
end
        
function generate_trees!(pt::PairTree, perm::Vector{Int})
    q = Int(floor(pt.n^(1/4)))
    row = col = Int(ceil(pt.n / q))
    for (x, y) in enumerate(perm)
        add_NDP(pt.tree1, [x-1, y], 1)
        add_NDP(pt.tree2, [x-1, y], sum_box_NDP(pt.tree1, [(0, x-1), (0, y)]))

        v = Int(ceil(x / q)) - 1

        h = Int(ceil((y + 1) / q)) - 1
        
        if !haskey(pt.Vert, v)
            pt.Vert[v] = Set{Tuple{Int,Int}}()
        end
        # println("v in gen trees = ", v)
        if !haskey(pt.Hori, h)
            pt.Hori[h] = Set{Tuple{Int,Int}}()
        end
        
        push!(pt.Vert[v], (x-1, y))
        push!(pt.Hori[h], (x-1, y))
    end
    
    for (x, y) in enumerate(perm)
        for s in 0:(row-1)
            if !haskey(pt.Vert_Trees, s)
                pt.Vert_Trees[s] = NDProductTree(pt.n, 2)
            end
            if !haskey(pt.Hori_Trees, s)
                pt.Hori_Trees[s] = NDProductTree(pt.n, 2)
            end
            
            add_NDP(pt.Vert_Trees[s], [x-1, y], 
                   sum_box_NDP(pt.tree1, [(0, min(x-1, (s + 1) * q)), (0, y)]))
    
            add_NDP(pt.Hori_Trees[s], [x-1, y], 
                   sum_box_NDP(pt.tree1, [(0, x-1), (0, min(y, (s + 1) * q))]))
        end
    end
end

function count_asc_pairs(pt::PairTree, perm::Vector{Int}, bounds)
    q = Int(floor(pt.n^(1/4)))
    
    a = Int(floor((bounds[1][1]) / q))
    b = Int(floor((bounds[1][2] - 1) / q))
    c = Int(floor((bounds[2][1]) / q))
    d = Int(floor((bounds[2][2] - 1)/ q))

    M = Set{Tuple{Int,Int}}()
    
    for i in [a, b]
        if haskey(pt.Vert, i)
            for v in pt.Vert[i]
                # println("v = ", v)
                if (bounds[1][1] <= v[1] < bounds[1][2]) && (bounds[2][1] <= v[2] < bounds[2][2])
                    push!(M, v)
                end
            end
        end
    end
    for j in [c, d]
        if haskey(pt.Hori, j)
            for h in pt.Hori[j]
                if (bounds[2][1] <= h[2] < bounds[2][2]) && (bounds[1][1] <= h[1] < bounds[1][2])
                    push!(M, h)
                end
            end
        end
    end

    asc_pairs = 0
    
    # ---------------------------------------------
    # Ascending pairs that end in M.
    # ---------------------------------------------

    for point in M
        asc_pairs += sum_box_NDP(pt.tree1, [(bounds[1][1], point[1]), (bounds[2][1], point[2])])
    end

    # ---------------------------------------------------
    # Ascending pairs that start in M and end in R_in.
    # ---------------------------------------------------

    for point in M
        if (point[1] < (a + 1) * q && point[2] < (c + 1) * q)
            asc_pairs += sum_box_NDP(pt.tree1, [((a + 1) * q, b * q), ((c + 1) * q, d * q)])
        elseif (point[1] < (a + 1) * q && (c + 1) * q <= point[2] < d * q)
            asc_pairs += sum_box_NDP(pt.tree1, [((a + 1) * q, b * q), (point[2] + 1, d * q)])
        elseif (point[2] < (c + 1) * q && (a + 1) * q <= point[1] < b * q)
            asc_pairs += sum_box_NDP(pt.tree1, [(point[1] + 1, b * q), ((c + 1) * q, d * q)])
        end
    end

    # --------------------------------
    # Ascending pairs in R_in.
    # --------------------------------

    tree2_sum = sum_box_NDP(pt.tree2, [((a + 1) * q, b * q), ((c + 1) * q, d * q)])
    vert_sum = sum_box_NDP(pt.Vert_Trees[a], [((a + 1) * q, b * q), ((c + 1) * q, d * q)])
    hori_sum = sum_box_NDP(pt.Hori_Trees[c], [((a + 1) * q, b * q), ((c + 1) * q, d * q)]) 
    corner_product = sum_box_NDP(pt.tree1, [((a + 1) * q, b * q), ((c + 1) * q, d * q)]) * 
                    sum_box_NDP(pt.tree1, [(0, (a + 1) * q), (0, (c + 1) * q)])

    asc_pairs += tree2_sum - vert_sum - hori_sum + corner_product
    return asc_pairs
end

# -------------------------------------------------------------------------------------------------
# Counts descending pairs by subtracting number of ascending pairs from total number of pairs.
# -------------------------------------------------------------------------------------------------

function count_desc_pairs(pt::PairTree, perm::Vector{Int}, bounds)
    asc = count_asc_pairs(pt, perm, bounds)
    total = binomial(sum_box_NDP(pt.tree1, bounds), 2)
    desc = total - asc
    return desc
end

"""
========================================================================================
Counts marked occurences of 4321(5) in a permutation. [Beniamini, Lavee]
========================================================================================
"""

function marked_count_43215(perm::Vector{Int})
    n = length(perm)
    chunk = Int(floor(n^(1/4)))
    permT = invperm(perm)
    T_out = NDProductTree(n, 2)

    # changing perm and permT to index 0
    perm = perm .-1
    permT = permT .-1

    #-------------------------------
    # 4,5 not in same row
    #-------------------------------

    for row in 0:chunk:n-1
        T1 = NDProductTree(n, 2)
        T2 = NDProductTree(n, 2)
        T3 = NDProductTree(n, 2)
        T4 = NDProductTree(n, 2)

        for i in 0:n-1
            if perm[i+1] < row
                add_NDP(T1, [i, perm[i+1]], 1)
                add_NDP(T2, [i, perm[i+1]], sum_box_NDP(T1, [(0, i), (perm[i+1], row)]))
                add_NDP(T3, [i, perm[i+1]], sum_box_NDP(T2, [(0, i), (perm[i+1], row)]))
                add_NDP(T4, [i, perm[i+1]], sum_box_NDP(T3, [(0, i), (perm[i+1], row)]))
            elseif row <= perm[i+1] < row + chunk
                add_NDP(T_out, [i, perm[i+1]], sum_box_NDP(T4, [(0, i), (0, row)]))
            end
        end
    end

    # -------------------------------------------------------------            
    # if 1,5 not in same col, but 4,5 in same row
    #--------------------------------------------------------------
   
    for row in chunk:chunk:n-1
        TT1 = NDProductTree(n, 2)
        TT2 = NDProductTree(n, 2)
        TT3 = NDProductTree(n, 2)
        TT4 = NDProductTree(n, 2)
        
        for (x, y) in enumerate(permT)
            if y < row
                add_NDP(TT1, [x-1, y], 1)
                add_NDP(TT2, [x-1, y], sum_box_NDP(TT1, [(0, x-1), (y, row)]))
                add_NDP(TT3, [x-1, y], sum_box_NDP(TT2, [(0, x-1), (y, row)]))
                add_NDP(TT4, [x-1, y], sum_box_NDP(TT3, [(0, x-1), (y, row)]))
            elseif row <= y < row + chunk
                col = (x-1) - (x - 1) % chunk
                add_NDP(T_out, [y, x-1], sum_box_NDP(TT4, [(col, x-1), (0, row)]))
            end
        end
    end

    # ------------------------------------------------
    # if 1,5 in same col and 4,5 in same row 
    # ------------------------------------------------

    PT = PairTree(n)
    generate_trees!(PT, perm)
    
    for (x5, y5) in enumerate(perm)
        col = (x5-1) - (x5 - 1) % chunk
        row = y5 - y5  % chunk
        
        for x1 in col:x5-1
            y1 = perm[x1+1]
            for y4 in row:y5-1
                x4 = permT[y4+1]
                if x4 < x1 && y4 > y1
                    add_NDP(T_out, [x5-1, y5], count_desc_pairs(PT, perm, [(x4+1, x1), (y1+1, y4)]))
                end
            end
        end
    end
                
    return sum_box_NDP(T_out, [(0, n), (0, n)])
end

rev1(p::Vector{Int}) = reverse(p)

rev2(p::Vector{Int}) = [(length(p)) - (x-1) for x in p]