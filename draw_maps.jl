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
                        -ts 800 400
                        -ot Byte
                        mapping/MAMMALS/MAMMALS.shp
                        $(fname)
                    `
            run(query)
            if iszero(filesize(fname))
                rm(fname)
            end
        catch e
            @info "nope"
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
        if lat + stride(msk,2) <= lc.bottom
            msk[lon,lat] = nothing
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

ric = convert(Float64, ric)
plot(ric)

# LCDB / SCBD
patches = findall(!isnothing, msk.grid)

# Them chonky bois are sparse
Y_host = spzeros(Int64, length(patches), length(hosts))
Y_virus_clover = spzeros(Int64, length(patches), length(viruses))
Y_clover = spzeros(Int64, length(patches), links(CLOVER))
Y_virus_trefle = spzeros(Int64, length(patches), length(viruses))
Y_trefle = spzeros(Int64, length(patches), links(TREFLE))

int_trefle = interactions(TREFLE)
int_clover = interactions(CLOVER)

for (i, tax) in enumerate(keys(ranges))
    @info tax
    istax = isequal(tax)
    sp_occ = findall(!isnothing, ranges[tax].grid)
    Y_host[indexin(sp_occ, patches), i] .= 1
    vir_pos_clover = indexin([x.from for x in filter(t -> istax(t.to), int_clover)], viruses)
    vir_pos_trefle = indexin([x.from for x in filter(t -> istax(t.to), int_trefle)], viruses)
    Y_clover[indexin(sp_occ, patches), findall(t -> istax(t.to), int_clover)] .= 1
    Y_virus_clover[indexin(sp_occ, patches), vir_pos] .= 1
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

LCBD(Yhost)
LCBD(Yvirus)
LCBD(Yint)

lcbd_host = similar(ranges[first(hosts)])
lcbd_virus = similar(ranges[first(hosts)])
lcbd_clover = similar(ranges[first(hosts)])

lcbd_host.grid[patches] = LCBD(Yhost)[1]
lcbd_virus.grid[patches] = LCBD(Yvirus)[1]
lcbd_clover.grid[patches] = LCBD(Yclover)[1]

p1 = heatmap(lcbd_host, frame=:box, legend=false)
p2 = heatmap(lcbd_virus, frame=:box, legend=false)
p3 = heatmap(lcbd_clover, frame=:box, legend=false)
for p in [p1, p2]
    yaxis!(p, (-60, 90), "Latitude")
    xaxis!(p, (-180, 180), "Longitude")
end
title!(p1, "Hosts")
title!(p2, "Viruses")

plot(p1, p2, layout=(2, 1), size=(400, 600))
savefig("lcbd_species.png")

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