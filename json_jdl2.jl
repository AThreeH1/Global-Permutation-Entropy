#!/usr/bin/env julia
###
# Convert the profile dictionaries from JSON to JLD2 format for faster loading in Julia.
###
using JSON, JLD2

# Hard-coded paths
JSON_PATH = "./profile_dicts/profile_6.json"
JLD2_PATH = "./profile_dicts/profile_6.jld2"

function reconstruct(x)
    # Dict marker?
    if x isa Dict{String,Any}
        if haskey(x, "__fraction__")
            n, d = x["__fraction__"]
            return n // d
        elseif haskey(x, "__tuple__")
            return tuple( reconstruct.(x["__tuple__"])... )
        elseif haskey(x, "__set__")
            return Set(reconstruct.(x["__set__"]))
        elseif haskey(x, "__dict__")
            d = Dict{Any,Any}()
            for pair in x["__dict__"]
                k = reconstruct(pair[1])
                v = reconstruct(pair[2])
                d[k] = v
            end
            return d
        else
            # plain JSON object
            d = Dict{String,Any}()
            for (k,v) in x
                d[k] = reconstruct(v)
            end
            return d
        end

    # Array → recurse
    elseif x isa Vector
        return reconstruct.(x)

    # primitive
    else
        return x
    end
end

# ——— Main ———
function gen_jdl2()
    # 1) Load JSON
    raw = JSON.parsefile(JSON_PATH)

    # 2) Reconstruct into pure Julia
    julia_profiles = reconstruct(raw)

    # 3) Save to JLD2
    @save JLD2_PATH julia_profiles

    println("JLD2 written to ", JLD2_PATH)

end

gen_jdl2()
