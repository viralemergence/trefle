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

function emptyraster()
    rast = geotiff(SimpleSDMResponse, joinpath(@__DIR__, "mask.tif"))
    rast.grid[findall(!isnothing, rast.grid)] .= 0.0
    return rast
end

# Richness (hosts)
richness = emptyraster()
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

# Number of zoonotic hosts
zhc = emptyraster()
zht = emptyraster()

# Number of zoonotic viruses
zvc = emptyraster()
zvt = emptyraster()

# Number of viruses
avc = emptyraster()
avt = emptyraster()

# Number of interactions
inc = emptyraster()
int = emptyraster()

zhc_names = filter(s -> s in keys(ranges), species(zCLOVER; dims=2))
for sp in zhc_names
    zhc.grid[findall(!isnothing, ranges[sp].grid)] .+= 1.0
end
zht_names = filter(s -> s in keys(ranges), species(zTREFLE; dims=2))
for sp in zht_names
    zht.grid[findall(!isnothing, ranges[sp].grid)] .+= 1.0
end

for i in findall(!isnothing, richness)
    if richness[i] > 0.0
        # Hosts at this location
        pool = [k for (k, v) in ranges if !isnothing(v[i])]
        # CLOVER
        N = simplify(CLOVER[:, pool])
        zvc[i] = count(x -> x in zoo_clover, species(N; dims=1))
        avc[i] = EcologicalNetworks.richness(N; dims=1)
        inc[i] = links(N)
        # TREFLE
        N = simplify(TREFLE[:, pool])
        zvt[i] = count(x -> x in zoo_trefle, species(N; dims=1))
        avt[i] = EcologicalNetworks.richness(N; dims=1)
        int[i] = links(N)
    end
end

function sprinkle(layer::T) where {T<:SimpleSDMLayer}
    return (longitudes(layer), latitudes(layer), transpose(replace(layer.grid, nothing => NaN)))
end

# Plot using GeoMakie
begin
    _proj = "natearth2"
    _coast = true
    _pal_known = :linear_wcmr_100_45_c42_n256
    _pal_gained = :linear_worb_100_25_c53_n256
    _pal_divergence = Reverse(:roma)

    fig = Figure(resolution=(1100, 780))

    ga1 = GeoAxis(fig[1, 2]; dest="+proj=$(_proj)", coastlines=_coast, title="A", subtitle="Known zoonotic hosts\n\n", titlealign=:left, subtitlecolor=:gray20)
    pl1 = GeoMakie.surface!(ga1, sprinkle(zhc)...; shading=false, interpolate=false, colormap=_pal_known)

    ga2 = GeoAxis(fig[2, 2]; dest="+proj=$(_proj)", coastlines=_coast, title="B", subtitle="Predicted new zoonotic hosts\n\n", titlealign=:left, subtitlecolor=:gray20)
    pl2 = GeoMakie.surface!(ga2, sprinkle(zht - zhc)...; shading=false, interpolate=false, colormap=_pal_gained)

    _hotspots = rescale(zht - zhc, (0, 1)) - rescale(richness, (0, 1))
    ga6 = GeoAxis(fig[3, 2]; dest="+proj=$(_proj)", coastlines=_coast, title="C", subtitle="Hotspots of zoonotic hosts gain\n\n", titlealign=:left, subtitlecolor=:gray20)
    pl6 = GeoMakie.surface!(ga6, sprinkle(_hotspots)...; shading=false, interpolate=false, colormap=_pal_divergence, colorrange=(-0.5, 0.5))

    ga3 = GeoAxis(fig[1, 3]; dest="+proj=$(_proj)", coastlines=_coast, title="D", subtitle="Known interactions\n\n", titlealign=:left, subtitlecolor=:gray20)
    pl3 = GeoMakie.surface!(ga3, sprinkle(inc)...; shading=false, interpolate=false, colormap=_pal_known)

    ga4 = GeoAxis(fig[2, 3]; dest="+proj=$(_proj)", coastlines=_coast, title="E", subtitle="Predicted new interactions\n\n", titlealign=:left, subtitlecolor=:gray20)
    pl4 = GeoMakie.surface!(ga4, sprinkle(int - inc)...; shading=false, interpolate=false, colormap=_pal_gained)

    _hotspots = rescale(int - inc, (0, 1)) - rescale(richness, (0, 1))
    ga5 = GeoAxis(fig[3, 3]; dest="+proj=$(_proj)", coastlines=_coast, title="F", subtitle="Hotspots of interaction gain\n\n", titlealign=:left, subtitlecolor=:gray20)
    pl5 = GeoMakie.surface!(ga5, sprinkle(_hotspots)...; shading=false, interpolate=false, colormap=_pal_divergence, colorrange=(-0.5, 0.5))

    cb1 = Colorbar(fig[1, 1], pl1; height=Relative(0.45))
    cb2 = Colorbar(fig[2, 1], pl2; height=Relative(0.45))
    cb6 = Colorbar(fig[3, 1], pl6; height=Relative(0.45))
    cb3 = Colorbar(fig[1, 4], pl3; height=Relative(0.45))
    cb4 = Colorbar(fig[2, 4], pl4; height=Relative(0.45))
    cb5 = Colorbar(fig[3, 4], pl5; height=Relative(0.45))

    datalims!(ga1)
    datalims!(ga2)
    datalims!(ga3)
    datalims!(ga4)
    datalims!(ga5)
    datalims!(ga6)

    fig
end
save(joinpath(@__DIR__, "..", "figures", "richness-pre-post.png"), fig, px_per_unit=2)

# LCDB / SCBD
@info "Finding the occupied patches"
patches = findall(!isnothing, richness.grid)

# Them chonky bois are not sparse anymore because beluga is our strong, robust son
@info "Allocating the arrays for LCBD"
Y_host = zeros(Int64, length(patches), length(hosts))
Y_virus_clover = zeros(Int64, length(patches), length(viruses))
Y_clover = zeros(Int64, length(patches), links(CLOVER))
Y_virus_trefle = zeros(Int64, length(patches), length(viruses))
Y_trefle = zeros(Int64, length(patches), links(TREFLE))

@info "Collecting interactions"
int_trefle = EcologicalNetworks.interactions(TREFLE)
int_clover = EcologicalNetworks.interactions(CLOVER)

@info "Preparing a list of species"
sp = collect(keys(ranges))

@info "Filling the LCBD arrays"
for i in 1:length(sp)
    tax = sp[i]
    @info tax
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
lcbd_host = emptyraster()
lcbd_virus_clover = emptyraster()
lcbd_virus_trefle = emptyraster()
lcbd_clover = emptyraster()
lcbd_trefle = emptyraster()

lcbd_host.grid[patches] = LCBD(Y_host)[1]
lcbd_virus_clover.grid[patches] = LCBD(Y_virus_clover)[1]
lcbd_virus_trefle.grid[patches] = LCBD(Y_virus_trefle)[1]
lcbd_clover.grid[patches] = LCBD(Y_clover)[1]
lcbd_trefle.grid[patches] = LCBD(Y_trefle)[1]

@info "Plotting"
begin
    _proj = "natearth2"
    _coast = true
    _pal_known = :linear_wcmr_100_45_c42_n256
    _pal_gained = :linear_worb_100_25_c53_n256
    _pal_divergence = Reverse(:roma)

    fig = Figure(resolution=(1100, 520))

    ga1 = GeoAxis(fig[1, 2]; dest="+proj=$(_proj)", coastlines=_coast, title="A", subtitle="Network uniqueness pre-imputation\n\n", titlealign=:left, subtitlecolor=:gray20)
    pl1 = GeoMakie.surface!(ga1, sprinkle(rescale(lcbd_clover, (0, 1)))...; shading=false, interpolate=false, colormap=_pal_known)

    ga2 = GeoAxis(fig[1, 3]; dest="+proj=$(_proj)", coastlines=_coast, title="B", subtitle="Network uniqueness post-imputation\n\n", titlealign=:left, subtitlecolor=:gray20)
    pl2 = GeoMakie.surface!(ga2, sprinkle(rescale(lcbd_trefle, (0, 1)))...; shading=false, interpolate=false, colormap=_pal_known)

    _hotspots = rescale(lcbd_trefle, (0, 1)) - rescale(lcbd_clover, (0, 1))
    ga3 = GeoAxis(fig[2, 2]; dest="+proj=$(_proj)", coastlines=_coast, title="C", subtitle="Hotspots of uniqueness gain\n\n", titlealign=:left, subtitlecolor=:gray20)
    pl3 = GeoMakie.surface!(ga3, sprinkle(_hotspots)...; shading=false, interpolate=false, colormap=_pal_divergence, colorrange=(-0.3, 0.3))

    ga4 = GeoAxis(fig[2, 3]; dest="+proj=$(_proj)", coastlines=_coast, title="D", subtitle="Host community uniqueness\n\n", titlealign=:left, subtitlecolor=:gray20)
    pl4 = GeoMakie.surface!(ga4, sprinkle(rescale(lcbd_host, (0, 1)))...; shading=false, interpolate=false, colormap=_pal_known)

    cb1 = Colorbar(fig[1, 1], pl1; height=Relative(0.45))
    cb2 = Colorbar(fig[1, 4], pl2; height=Relative(0.45))
    cb3 = Colorbar(fig[2, 1], pl3; height=Relative(0.45))
    cb4 = Colorbar(fig[2, 4], pl4; height=Relative(0.45))

    datalims!(ga1)
    datalims!(ga2)
    datalims!(ga3)
    datalims!(ga4)

    fig
end
save(joinpath(@__DIR__, "..", "figures", "lcbd-pre-post.png"), fig, px_per_unit=2)