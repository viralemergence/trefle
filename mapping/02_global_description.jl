using SimpleSDMLayers
using Plots

using DataFrames
import CSV

using Statistics

# Load the rasters
ranges_path = joinpath(@__DIR__, "cleaned_rasters")
ranges_rasters = filter(endswith(".tif"), readdir(ranges_path))

# Transform the file path to a dictionary key
function file_to_key(f)
    return replace(first(split(f, ".")), "_" => " ")
end

# Layer of rasters
ranges = Dict([file_to_key(f) => geotiff(SimpleSDMPredictor, joinpath(ranges_path, f)) for f in ranges_rasters])

# Richness (hosts)
richness = geotiff(SimpleSDMResponse, joinpath(@__DIR__, "mask.tif"))
richness.grid[findall(!isnothing, richness.grid)] .= 0.0
for (k, v) in ranges
    richness.grid[findall(!isnothing, v.grid)] .+= 1.0
end

plot(richness)

# Viral richness
trefle = DataFrame(CSV.File(joinpath(@__DIR__, "..", "artifacts", "trefle.csv")))
clover = DataFrame(CSV.File(joinpath(@__DIR__, "..", "data", "clover.csv")))
hosts = unique(trefle.host)
viruses = unique(trefle.virus)

# Convert into networks
A = zeros(Bool, (length(viruses), length(hosts)))
TREFLE = BipartiteNetwork(A, viruses, hosts)
CLOVER = BipartiteNetwork(A, viruses, hosts)
for v in eachrow(trefle)
    TREFLE[v.virus, v.host] = true
end
for v in eachrow(clover)
    CLOVER[v.Virus, v.Host] = true
end

# Get the zoonotic components
zoo_clover = CLOVER[:, "Homo sapiens"]
zoo_trefle = TREFLE[:, "Homo sapiens"]

zCLOVER = simplify(CLOVER[collect(zoo_clover), :])
zTREFLE = simplify(TREFLE[collect(zoo_trefle), :])

Y = zeros(Int64, length(richness), length(ranges))
patches = findall(!isnothing, richness)
for (i, r) in enumerate(ranges)
    idx = indexin(findall(!isnothing, r.second), patches)
    Y[idx, i] .= 1
end

function hellinger(Y::Matrix{T}) where {T<:Number}
    yi = sum(Y; dims=2)
    return sqrt.(Y ./ yi)
end

function LCBD(Y)
    S = (Y .- mean(Y; dims=1)) .^ 2.0
    SStotal = sum(S)
    BDtotal = SStotal / (size(Y, 1) - 1)
    SSj = sum(S; dims=1)
    SCBDj = SSj ./ SStotal
    SSi = sum(S; dims=2)
    LCBDi = SSi ./ SStotal
    return LCBDi, SCBDj, BDtotal
end

host_lcbd = geotiff(SimpleSDMResponse, joinpath(@__DIR__, "mask.tif"))
host_lcbd.grid[findall(!isnothing, host_lcbd.grid)] .= LCBD(Y)[1]
plot(host_lcbd)