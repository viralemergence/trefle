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
imputation_rank = parse(Int64, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))

# Get them matrices
A = Array(clover.edges)

# Pools
filled = findall(A)
empty = findall(.!(A))

# Positions to sample
n_samples = Threads.nthreads()*20

dfs = [DataFrame(value = Bool[], position = [], model = Symbol[], initial = Float64[], updated = Float64[]) for i in 1:Threads.nthreads()]

QUER = [Float64.(A) for i in 1:Threads.nthreads()]
TMPLn1 = [n1(Q) for Q in QUER]
TMPLn2 = [n2(QUER[i]) for i in 1:Threads.nthreads()]
TMPLn3 = [n3(QUER[i]) for i in 1:Threads.nthreads()]

empties = rand(empty, n_samples)
filleds = rand(filled, n_samples)

Threads.@threads for i in 1:n_samples
    Qu = QUER[Threads.threadid()]
    try
        thisone = impute(Qu, TMPLn1[Threads.threadid()], empties[i]; rank=imputation_rank)
        push!(dfs[Threads.threadid()], (false, i, :connectance, thisone.initial, thisone.updated))
        thisone = impute(Qu, TMPLn2[Threads.threadid()], empties[i]; rank=imputation_rank)
        push!(dfs[Threads.threadid()], (false, i, :degree, thisone.initial, thisone.updated))
        thisone = impute(QUER[Threads.threadid()], TMPLn3[Threads.threadid()], empties[i]; rank=imputation_rank)
        push!(dfs[Threads.threadid()], (false, i, :hybrid, thisone.initial, thisone.updated))

        thisone = impute(Qu, TMPLn1[Threads.threadid()], filleds[i]; rank=imputation_rank)
        push!(dfs[Threads.threadid()], (true, i, :connectance, thisone.initial, thisone.updated))
        thisone = impute(Qu, TMPLn2[Threads.threadid()], filleds[i]; rank=imputation_rank)
        push!(dfs[Threads.threadid()], (true, i, :degree, thisone.initial, thisone.updated))
        thisone = impute(Qu, TMPLn3[Threads.threadid()], filleds[i]; rank=imputation_rank)
        push!(dfs[Threads.threadid()], (true, i, :hybrid, thisone.initial, thisone.updated))
    catch
        continue
    end
    CSV.write("svd-tuning-thread-$(Threads.threadid())-rank-$(imputation_rank).csv", dfs[Threads.threadid()])
end
