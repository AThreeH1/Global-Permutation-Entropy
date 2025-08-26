using Pkg
using PyCall

packages = [
    "Combinatorics",
    "Plots",
    "Symbolics",
    "BenchmarkTools",
    "JLD2",
    "PyCall",
    "JSON",
]

for pkg in packages
    println("Installing or updating package: $pkg")
    Pkg.add(pkg)
end

println("All packages are installed.")

# Get the directory containing this script (where pkl_json.py is located)
script_dir = dirname(abspath(@__FILE__))

# Add directory to Python's sys.path
sys = pyimport("sys")
if script_dir ∉ sys."path"
    pushfirst!(PyVector(sys."path"), script_dir)
end

# Import the pkl_json module
pkl_json = pyimport("pkl_json")

# Check for the file using absolute path
filepath = joinpath(script_dir, "profile_dicts", "profile_6.jld2")
println("Looking for file at: ", filepath)
println("File exists: ", isfile(filepath))

if !isfile(filepath)
    println("File not found. Running generators...")
    pkl_json.gen_json()
    include("json_jdl2.jl")
    gen_jdl2()
    println("After generation, file exists: ", isfile(filepath))
else
    println("File exists. Skipping generation.")
end

println("To run the script in parallel, use the following command: JULIA_NUM_THREADS=32 julia PATH_TO_FILE.jl")

