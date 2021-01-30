using DataFrames
import CSV
using EcologicalNetworks
using StatsPlots
using LinearAlgebra
using Statistics

# Functions
include("lib.jl")

# Read clover
include("read_clover.jl")

# Get array parameters
virus_to_predict = parse(Int64, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
imputation_rank = 12

# Get them matrices
A = Array(clover.edges)

# Pools
filled = findall(A)
empty = findall(.!(A))

dfs = [DataFrame(value = Bool[], host = String[], virus = String[], host_pos = Int64[], virus_pos = Int64[], initial = Float64[], updated = Float64[]) for i in 1:Threads.nthreads()]

QUER = [Float64.(A) for i in 1:Threads.nthreads()]
TMPL = [n1(Q) for Q in QUER]

@time Threads.@threads for i in 1:size(clover, 2)
    Qu = QUER[Threads.threadid()]
    try
        pos = CartesianIndex(virus_to_predict, i)
        result = impute(Qu, TMPL[Threads.threadid()], pos; rank=12)
        push!(dfs[Threads.threadid()],
            (false, hosts[i], viruses[virus_to_predict], i, virus_to_predict, result.initial, result.updated)
        )
    catch
        @warn "Unable to do prediction for $(viruses[virus_to_predict]) and host $(hosts[i])"
        continue
    end
    CSV.write("svd-prediction-thread-$(Threads.threadid())-virus-$(virus_to_predict).csv", dfs[Threads.threadid()])
end