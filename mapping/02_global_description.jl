using SimpleSDMLayers
using EcologicalNetworks

using DataFrames
import CSV

using Statistics

# Plots with spatial wrapping
using GeoMakie, CairoMakie

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

# Number of zoonotic hosts/viruses
zhc = geotiff(SimpleSDMResponse, joinpath(@__DIR__, "mask.tif"))
zhc.grid[findall(!isnothing, zhc.grid)] .= 0.0
zht = geotiff(SimpleSDMResponse, joinpath(@__DIR__, "mask.tif"))
zht.grid[findall(!isnothing, zht.grid)] .= 0.0

zhc_names = filter(s -> s in keys(ranges), species(zCLOVER; dims=2))
for sp in zhc_names
    zhc.grid[findall(!isnothing, ranges[sp].grid)] .+= 1.0
end
zht_names = filter(s -> s in keys(ranges), species(zTREFLE; dims=2))
for sp in zht_names
    zht.grid[findall(!isnothing, ranges[sp].grid)] .+= 1.0
end

zvc = geotiff(SimpleSDMResponse, joinpath(@__DIR__, "mask.tif"))
zvc.grid[findall(!isnothing, zvc.grid)] .= 0.0
zvt = geotiff(SimpleSDMResponse, joinpath(@__DIR__, "mask.tif"))
zvt.grid[findall(!isnothing, zvt.grid)] .= 0.0
for i in findall(!isnothing, richness)
    if richness[i] > 0.0
        pool = [k for (k, v) in ranges if !isnothing(v[i])]
        # CLOVER
        N = simplify(CLOVER[:, pool])
        v = filter(x -> x in zoo_clover, species(N; dims=1))
        zvc[i] = length(v)
        # TREFLE
        N = simplify(TREFLE[:, pool])
        v = filter(x -> x in zoo_trefle, species(N; dims=1))
        zvt[i] = length(v)
    end
end

# Plot using GeoMakie
_proj = "natearth"
fig = Figure(resolution=(900, 700))
ga1 = GeoAxis(fig[1, 2]; dest="+proj=$(_proj)", coastlines=true, title="Zoonotic hosts pre-imputation\n\n")
pl1 = GeoMakie.surface!(ga1, longitudes(richness), latitudes(richness), transpose(replace(zhc.grid, nothing => NaN)); shading=false, interpolate=false, colormap=:linear_wcmr_100_45_c42_n256)
ga2 = GeoAxis(fig[2, 2]; dest="+proj=$(_proj)", coastlines=true, title="Zoonotic hosts gained post-imputation\n\n")
pl2 = GeoMakie.surface!(ga2, longitudes(richness), latitudes(richness), transpose(replace((zht - zhc).grid, nothing => NaN)); shading=false, interpolate=false, colormap=:linear_worb_100_25_c53_n256)
cb1 = Colorbar(fig[1, 1], pl1; height=Relative(0.55))
cb2 = Colorbar(fig[2, 3], pl2; height=Relative(0.55))
datalims!(ga1)
datalims!(ga2)
fig
save("zoo-map-draft.png", fig, px_per_unit=2)

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
