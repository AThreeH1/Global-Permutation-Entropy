# Global-Permutation-Entropy
Package to calculate the complexity measure Global Permutation Entropy
introduced in (REF).

# Setup
Run the setup. This installs the necessary packages and generates the necessary files in profile_dicts/.
 ```bash
    $ julia setup.jl
 ```


# Example
Import GPE.jl and use global_permutation_entropy to calculate entropy.
For example:
```julia
    include("GPE.jl")
    random_data = randn(20)
    entropy = global_permutation_entropy(random_data, order=3, window_size=10)
    println(entropy) # Length 11, since there are 11 windows of size 10. Should give quite high entropy values for this random data.

    strictly_increasing_data = collect(1:20)
    println(strictly_increasing_data)
    entropy = global_permutation_entropy(strictly_increasing_data, order=3, window_size=10)
    println(entropy) # Length 11, should give entropy 0. since there is only one permutation pattern appearing.
```

# Tests
To run the tests
```bash
    $ julia test_profiles.jl
```

# Parallel
To use parallel execution, one needs to run julia with 
the JULIA_NUM_THREADS=32 environment variable set.
For example, to use 4 threads in the REPL:
```bash
    $ JULIA_NUM_THREADS=4 julia
```
or to run the tests:
```bash
    $ JULIA_NUM_THREADS=4 julia test_profiles.jl
```


# Implementation Details:

The profiles are implemented in Julia, leveraging its built-in support for parallel computing allowing for efficient execution on multi-core systems.
The core data structure is a dictionary, where each key is labeled as a specific algorithm
class (e.g., corner trees, tree double posets, count_NE_patterntree, etc), and each value
stores the associated computational data—such as tree structures with edge constraints or permutation patterns—required for enumeration. The idea is to enumerate through this dictionary and have conditional statements that feeds the value to the appropriate algorithms.

We employ two levels of parallelism, External over data batches and one Internal over algorithm classes, depending on the profile order:


1. External Parallelization (Profiles 2 and 3):
For smaller profiles, the number of structures to be counted is relatively
small, so internal threading over algorithm classes yields limited speedup.
Instead, we parallelize across data batches. Once the user specifies a signal,
window size, and stride, the data is segmented into
(overlapping or disjoint) windows (batches).
Each thread computes the profile independently on a
separate batch. This method efficiently utilizes multiple cores without
incurring significant overhead from task granularity.

2. Internal Parallelization (Profiles 4, 5, and 6):  
For larger profiles, the number of distinct structures 213 pattern trees
in the 6-profile increases substantially. Here, we group structures by their
algorithm type - corner trees, pattern trees, or tree double posets - and apply
fine-grained parallelism within each group. For instance, in the 5-profile,
all n pattern trees are evaluated in parallel across threads. This internal parallelization yields significant
performance gains, especially since these computations dominate the runtime. The results from each threaded task are aggregated after computation.
No external parallelization is used here.