# Global-Permutation-Entropy
Complexity measure for ordinal patterns

To calculate GPE:
 Step 1:
 Run setup.jl
 Step 2:
 Import GPE.jl and use global_permutation_entropy to calculate entropy. 

Implementation Details:

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




