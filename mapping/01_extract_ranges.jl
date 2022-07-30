using DataFrames
import CSV

using SimpleSDMLayers
using EcologicalNetworks

import GDAL
using ArchGDAL

# Load the trefle / clover edgelists
trefle = DataFrame(CSV.File(joinpath(@__DIR__, "..", "artifacts", "trefle.csv")))
clover = DataFrame(CSV.File(joinpath(@__DIR__, "..", "data", "clover.csv")))
hosts = unique(trefle.host)
viruses = unique(trefle.virus)

# Prepare a directory to store the rasters
raster_path = joinpath(@__DIR__, "..", "mapping", "rasters")
isdir(raster_path) || mkdir(raster_path)

cleaned_raster_path = joinpath(@__DIR__, "..", "mapping", "cleaned_rasters")
isdir(cleaned_raster_path) || mkdir(cleaned_raster_path)

# Prepare the IUCN data
iucn_path = get(ENV, "IUCN_PATH", nothing)
if isnothing(iucn_path)
    @error "You need an environmental variable `IUCN_PATH` pointing to the rangemaps"
end

# Check that the MAMMALS data are present
if ~isdir(joinpath(iucn_path, "MAMMALS"))
    @error "You need the unzipped `MAMMALS` folder in the IUCN_PATH"
end

# Download landcover data to remove the pixels in full open water
open_water = SimpleSDMPredictor(EarthEnv, LandCover, 12)

# Read the layers
for host in hosts[1:5]
    @info "Rasterizing $(host)"
    host_path = joinpath(raster_path, replace(host, " " => "_") * ".tif")
    if ~isfile(host_path)
        try
            rasterize_query = `gdal_rasterize
                -l "MAMMALS"
                -where "binomial = '$(host)'"
                -a presence
                -te $(open_water.left) $(open_water.bottom) $(open_water.right) $(open_water.top)
                -ts 200 100
                -ot Byte
                -add
                $(joinpath(iucn_path, "MAMMALS", "MAMMALS.shp"))
                $(host_path)`
            run(rasterize_query)
        catch err
            print(err)
        end
    end
end

# Function to coerce the data from a layer to the other
function nukepix(reference::TR, target::TT, f) where {TR<:SimpleSDMLayer,TT<:SimpleSDMLayer}
    coerced = SimpleSDMResponse(zeros(Bool, size(target)), target)
    for k in keys(coerced)
        coerced[k] = f(clip(reference, k .- stride(reference), k .+ stride(reference)))
    end
    return coerced
end

# Make a mask from open water pixels
first_raster = first(readdir(raster_path; join=true))
range_ref = geotiff(SimpleSDMPredictor, first_raster)
coerced = nukepix(open_water, range_ref, g -> ~all(g.grid .== 100))
replace!(coerced.grid, false => nothing)

# Cleanup the rasterized data
for layer in readdir(raster_path)
    this_raster = joinpath(raster_path, layer)
    @info "Cleaning up $(layer)"
    if ~iszero(filesize(this_raster))
        l = geotiff(SimpleSDMPredictor, this_raster)
        l.grid[findall(v -> ~isnothing(v), l.grid)] .= one(eltype(l))
        clean_layer = convert(Float64, mask(coerced, l))
        geotiff(joinpath(cleaned_raster_path, layer), clean_layer)
    end
end
