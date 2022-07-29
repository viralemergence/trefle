using SimpleSDMLayers
using Plots
using DataFrames
import CSV
using EcologicalNetworks
using Statistics
using StatsBase
using GLM
using SparseArrays
import GDAL
using ArchGDAL

# List of species in trefle
trefle = DataFrame(CSV.File(joinpath(@__DIR__, "artifacts", "trefle.csv")))
clover = DataFrame(CSV.File(joinpath(@__DIR__, "data", "clover.csv")))
hosts = unique(trefle.host)
viruses = unique(trefle.virus)

# Interactions
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

raster_path = joinpath(@__DIR__, "mapping", "rasters")
ispath(raster_path) || mkdir(raster_path)

# Get the predictors
ranges = Dict{String,SimpleSDMPredictor}()
for host in hosts
    @info host
    fname = joinpath(raster_path, replace(host, " " => "_") * ".tif")
    if ~isfile(fname)
        try
            query = `gdal_rasterize
                        -l "MAMMALS"
                        -where "binomial = '$(host)'"
                        -a presence
                        -ts 300 150
                        -ot Byte
                        mapping/MAMMALS/MAMMALS.shp
                        $(fname)
                    `
            run(query)
        catch e
            @info "nope"
        end
        if iszero(filesize(fname))
            rm(fname)
        end
    else
        if ~iszero(filesize(fname))
            mp = SimpleSDMLayers.geotiff(SimpleSDMPredictor, fname, 1)
            replace!(mp.grid, zero(eltype(mp)) => nothing)
            mp.grid[findall(v -> ~isnothing(v), mp.grid)] .= one(eltype(mp))
            ranges[host] = mp
            mp = nothing
        end
    end
    GC.gc()
end

# Make a mask to remove fully open water pixels
msk = convert(Bool, similar(ranges[first(keys(ranges))]))
msk.grid[:, :] .= nothing
lc = SimpleSDMPredictor(EarthEnv, LandCover, 12)
for lat in latitudes(msk)
    for lon in longitudes(msk)
        try
            tlc = clip(lc; left=lon - 1.001stride(msk, 1), right=lon + 1.001stride(msk, 1), bottom=lat - 1.001stride(msk, 2), top=lat + 1.001stride(msk, 2))
            if ~all(isequal(100), tlc.grid)
                msk[lon, lat] = true
            end
        catch e
        end
        if lat + stride(msk, 2) <= lc.bottom
            msk[lon, lat] = nothing
        end
        if lon + stride(msk, 1) >= lc.right
            msk[lon, lat] = nothing
        end
    end
end

maskedranges = Dict([s => mask(msk, ranges[s]) for s in keys(ranges)])

ric = similar(ranges[first(keys(maskedranges))])
for (k, v) in maskedranges
    for i in findall(!isnothing, v.grid)
        if isnothing(ric.grid[i])
            ric.grid[i] = 1
        else
            ric.grid[i] += 1
        end
    end
end

richness = convert(Float64, ric)
plot(richness, frame=:box)

# LCDB / SCBD
@info "Finding the occupied patches"
patches = findall(!isnothing, richness.grid)

# Them chonky bois are not sparse anymore because beluga is our strong, robust son
@info "Allocating the arrays for LCBD"
Y_host = spzeros(Int64, length(patches), length(hosts))
Y_virus_clover = spzeros(Int64, length(patches), length(viruses))
Y_clover = spzeros(Int64, length(patches), links(CLOVER))
Y_virus_trefle = spzeros(Int64, length(patches), length(viruses))
Y_trefle = spzeros(Int64, length(patches), links(TREFLE))

@info "Collecting interactions"
int_trefle = interactions(TREFLE)
int_clover = interactions(CLOVER)

@info "Preparing a list of species"
sp = collect(keys(ranges))

@info "Filling the LCBD arrays"
Threads.@threads for i in 1:length(sp)
    tax = sp[i]
    try
        istax = isequal(tax)
        sp_occ = findall(!isnothing, ranges[tax].grid)
        Y_host[indexin(sp_occ, patches), i] .= 1
        vir_pos_clover = indexin([x.from for x in filter(t -> istax(t.to), int_clover)], viruses)
        vir_pos_trefle = indexin([x.from for x in filter(t -> istax(t.to), int_trefle)], viruses)
        Y_clover[indexin(sp_occ, patches), findall(t -> istax(t.to), int_clover)] .= 1
        Y_trefle[indexin(sp_occ, patches), findall(t -> istax(t.to), int_trefle)] .= 1
        Y_virus_clover[indexin(sp_occ, patches), vir_pos_clover] .= 1
        Y_virus_trefle[indexin(sp_occ, patches), vir_pos_trefle] .= 1
    catch e
    end
end

@info "Declaring LCBD functions"
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

# Raw version
@info "Pre-allocating response layers"
lcbd_host = similar(ranges[first(hosts)])
lcbd_virus_clover = similar(ranges[first(hosts)])
lcbd_virus_trefle = similar(ranges[first(hosts)])
lcbd_clover = similar(ranges[first(hosts)])
lcbd_trefle = similar(ranges[first(hosts)])

lcbd_host.grid[patches] = LCBD(Y_host)[1]
lcbd_virus_clover.grid[patches] = LCBD(Y_virus_clover)[1]
lcbd_virus_trefle.grid[patches] = LCBD(Y_virus_trefle)[1]
lcbd_clover.grid[patches] = LCBD(Y_clover)[1]
lcbd_trefle.grid[patches] = LCBD(Y_trefle)[1]

@info "Saving response layers"
SimpleSDMLayers.ascii(lcbd_host, "lcbd_host.ascii")
SimpleSDMLayers.ascii(lcbd_virus_clover, "lcbd_virus_clover.ascii")
SimpleSDMLayers.ascii(lcbd_virus_trefle, "lcbd_virus_trefle.ascii")
SimpleSDMLayers.ascii(lcbd_clover, "lcbd_clover.ascii")
SimpleSDMLayers.ascii(lcbd_trefle, "lcbd_trefle.ascii")
SimpleSDMLayers.ascii(richness, "host_richness.ascii")

# Richness map
@info "Viral richness maps"
vr_clover = similar(richness)
vr_trefle = similar(richness)

vr_clover.grid[patches] = vec(sum(Y_virus_clover; dims=2))
vr_trefle.grid[patches] = vec(sum(Y_virus_trefle; dims=2))

SimpleSDMLayers.ascii(vr_clover, "viral_richness_clover.ascii")
SimpleSDMLayers.ascii(vr_trefle, "viral_richness_trefle.ascii")

CSV.write("lcbd_host.csv", filter(r -> !isnothing(r.values), DataFrame(lcbd_host)))
CSV.write("lcbd_virus_clover.csv", filter(r -> !isnothing(r.values), DataFrame(lcbd_virus_clover)))
CSV.write("lcbd_virus_trefle.csv", filter(r -> !isnothing(r.values), DataFrame(lcbd_virus_trefle)))
CSV.write("lcbd_clover.csv", filter(r -> !isnothing(r.values), DataFrame(lcbd_clover)))
CSV.write("lcbd_trefle.csv", filter(r -> !isnothing(r.values), DataFrame(lcbd_trefle)))
CSV.write("richness_hosts.csv", filter(r -> !isnothing(r.values), DataFrame(richness)))
CSV.write("richness_virus_clover.csv", filter(r -> !isnothing(r.values), DataFrame(vr_clover)))
CSV.write("richness_virus_trefle.csv", filter(r -> !isnothing(r.values), DataFrame(vr_trefle)))

#=
ols = lm(@formula(Y ~ X), data)

X = Float64.(filter(!isnothing, lcbd_betacov_hosts.grid))
y = Float64.(filter(!isnothing, icbd_betacov_hosts.grid))

regmodel = lm(permutedims(permutedims(X)), y)

hotspots = similar(icbd_betacov_hosts)
hotspots.grid[patches] .= residuals(regmodel)

heatmap(hotspots, dpi=400, frame=:box, clim=(-0.015, 0.015), c=:PuOr)
yaxis!((-90,90), "Latitude")
xaxis!((-180,180), "Longitude")
savefig("betacov-residuals.png")
=#