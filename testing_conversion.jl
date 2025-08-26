#!/usr/bin/env julia

using PyCall, JLD2

function main()
    # — Paths — 
    PKL_PATH  = "./parallel_implementation/profile_dicts/profile_4.pkl"
    JLD2_PATH = "./parallel_implementation/profile_dicts/profile_4.jld2"

    # — Load JLD2 data —
    println("Loading JLD2 from ", JLD2_PATH)
    @load JLD2_PATH julia_profiles

    # — Load pickle via PyCall —
    println("Loading pickle via PyCall from ", PKL_PATH)
    pickle   = pyimport("pickle")
    builtins = pyimport("builtins")
    py_open  = builtins.open

    f = py_open(PKL_PATH, "rb")
    py_data = pickle.load(f)

    # — Prepare Python built-in type references —
    Fraction = pyimport("fractions").Fraction
    py_int   = pybuiltin("int")
    py_float = pybuiltin("float")
    py_str   = pybuiltin("str")
    py_bool  = pybuiltin("bool")
    py_dict  = pybuiltin("dict")
    py_list  = pybuiltin("list")
    py_tuple = pybuiltin("tuple")
    py_set   = pybuiltin("set")

    # — Minimal conversion from PyObject containers to Julia, leaving unknown numeric PyObjects intact —
    function py_to_julia(x)
        if !(x isa PyObject)
            return x
        elseif PyCall.isinstance(x, Fraction)
            # Python fractions.Fraction → Julia Rational
            num = Int(x[:numerator])
            den = Int(x[:denominator])
            return num // den
        elseif PyCall.isinstance(x, py_int)
            return Int(x)
        elseif PyCall.isinstance(x, py_float)
            return Float64(x)
        elseif PyCall.isinstance(x, py_str)
            # Convert Python str to Julia String via pycall
            return pycall(py_str, String, x)
        elseif PyCall.isinstance(x, py_bool)
            return Bool(x)
        elseif PyCall.isinstance(x, py_dict)
            d = Dict{Any,Any}()
            for k in x.keys()
                d[py_to_julia(k)] = py_to_julia(x[k])
            end
            return d
        elseif PyCall.isinstance(x, py_list)
            return [ py_to_julia(el) for el in x ]
        elseif PyCall.isinstance(x, py_tuple)
            arr = collect(x)
            return tuple(py_to_julia.(arr)...)
        elseif PyCall.isinstance(x, py_set)
            return Set(py_to_julia(el) for el in x)
        else
            # Leave as PyObject if not recognized built-in numeric or container
            return x
        end
    end

    println("Converting pickle data to pure Julia (containers unwrapped)…")
    data_from_pkl = py_to_julia(py_data)

    # Ensure both are indexable collections with same length
    println("Comparing lengths…")
    len_pkl = try
        length(data_from_pkl)
    catch
        error("Top-level PKL object is not indexable with length. Type: $(typeof(data_from_pkl))")
    end
    len_jld = try
        length(julia_profiles)
    catch
        error("Top-level JLD2 object is not indexable with length. Type: $(typeof(julia_profiles))")
    end

    if len_pkl != len_jld
        error("Length mismatch: PKL has $len_pkl entries, JLD2 has $len_jld entries")
    end

    println("Scanning for first mismatch…")

    # Comparison function: parse PyObject via Python str() and compare to Julia value
    function equal_pyobj_vs_julia(x, y)
        # 1. Direct equality: covers many Julia-only cases or identical types
        try
            if x == y
                return true
            end
        catch
            # comparing PyObject with Julia type might error; skip to parsing logic
        end

        # 2. Handle dictionary comparison recursively
        if isa(x, Dict) && isa(y, Dict)
            # Check if keys are the same
            keys_x = Set(keys(x))
            keys_y = Set(keys(y))
            if keys_x != keys_y
                return false
            end
            
            # Check if all values are equal recursively
            for k in keys_x
                if !equal_pyobj_vs_julia(x[k], y[k])
                    return false
                end
            end
            return true
        end

        # 3. Handle array/vector comparison recursively
        if (isa(x, AbstractArray) && isa(y, AbstractArray)) || 
           (isa(x, Tuple) && isa(y, Tuple))
            if length(x) != length(y)
                return false
            end
            for i in 1:length(x)
                if !equal_pyobj_vs_julia(x[i], y[i])
                    return false
                end
            end
            return true
        end

        # Helper: get Python str(x) as Julia String, safely
        get_py_str = function(obj)
            # obj is PyObject; call Python str(obj), then convert to Julia String
            try
                return pycall(py_str, String, obj)
            catch
                return nothing
            end
        end

        # 4. If x is PyObject, get s = str(x) and attempt to parse as integer or rational
        if x isa PyObject
            s = get_py_str(x)
            if s !== nothing
                # Try integer pattern: optional leading +/-, digits only
                if occursin(r"^[+-]?\d+$", s)
                    try
                        xi = parse(Int, s)
                        # Compare xi to y
                        if y isa Rational{<:Integer}
                            if xi == y
                                return true
                            end
                        elseif y isa Integer
                            if xi == y
                                return true
                            end
                        elseif y isa AbstractFloat
                            if Float64(xi) == y
                                return true
                            end
                        end
                    catch
                        # parse failed unexpectedly
                    end
                end
                # Try rational pattern: "num/den"
                if occursin(r"^[+-]?\d+/\d+$", s)
                    parts = split(s, '/')
                    if length(parts) == 2
                        p_str = parts[1]
                        q_str = parts[2]
                        try
                            p = parse(Int, p_str)
                            q = parse(Int, q_str)
                            if q != 0
                                pr = p // q
                                if y isa Rational{<:Integer}
                                    if pr == y
                                        return true
                                    end
                                elseif y isa Integer
                                    if pr == y
                                        return true
                                    end
                                elseif y isa AbstractFloat
                                    if Float64(pr) == y
                                        return true
                                    end
                                end
                            end
                        catch
                            # parse failed; skip
                        end
                    end
                end
            end
        end

        # 5. If y is PyObject, do symmetric parsing
        if y isa PyObject
            s = get_py_str(y)
            if s !== nothing
                if occursin(r"^[+-]?\d+$", s)
                    try
                        yi = parse(Int, s)
                        if x isa Rational{<:Integer}
                            if yi == x
                                return true
                            end
                        elseif x isa Integer
                            if yi == x
                                return true
                            end
                        elseif x isa AbstractFloat
                            if Float64(yi) == x
                                return true
                            end
                        end
                    catch
                    end
                end
                if occursin(r"^[+-]?\d+/\d+$", s)
                    parts = split(s, '/')
                    if length(parts) == 2
                        p_str = parts[1]
                        q_str = parts[2]
                        try
                            p = parse(Int, p_str)
                            q = parse(Int, q_str)
                            if q != 0
                                pr = p // q
                                if x isa Rational{<:Integer}
                                    if pr == x
                                        return true
                                    end
                                elseif x isa Integer
                                    if pr == x
                                        return true
                                    end
                                elseif x isa AbstractFloat
                                    if Float64(pr) == x
                                        return true
                                    end
                                end
                            end
                        catch
                        end
                    end
                end
            end
        end

        return false
    end

    # Main loop: compare each element
    for i in 1:len_pkl
        a = data_from_pkl[i]
        b = julia_profiles[i]

        if equal_pyobj_vs_julia(a, b)
            continue
        end

        # Not equal: report mismatch
        println("\n  MISMATCH at index $i:")
        println("→ From PKL: ")
        show(stdout, MIME("text/plain"), a)
        println("\n   typeof PKL value: ", typeof(a))
        println("\n→ From JLD2:")
        show(stdout, MIME("text/plain"), b)
        println("\n   typeof JLD2 value: ", typeof(b))

        # If both Dicts, diff keys/values using same comparison
        if isa(a, Dict) && isa(b, Dict)
            Ka = Set(keys(a))
            Kb = Set(keys(b))
            println("\n Keys only in PKL: ", collect(setdiff(Ka, Kb)))
            println(" Keys only in JLD2: ", collect(setdiff(Kb, Ka)))

            common = intersect(collect(Ka), collect(Kb))
            for k in common
                va = a[k]
                vb = b[k]
                if !equal_pyobj_vs_julia(va, vb)
                    println("\n • Key = ", k)
                    println("    PKL val = ", va, "  typeof=", typeof(va))
                    println("    JLD2 val= ", vb, "  typeof=", typeof(vb))
                end
            end
        end

        error("Stopped at first mismatch.")
    end

    println("Perfect match: all entries align (allowing PyObject N equal to N//1 via string parsing)!")
end

# Run main()
main()
