# Load the rasters
raster_path = joinpath(@__DIR__, "cleaned_rasters")
rasters = filter(endswith(".tif"), readdir(raster_path))

# Transform the file path to a dictionary key
function file_to_key(f)
    return replace(first(split(f, ".")), "_" => " ")
end

# Layer of rasters
ranges = Dict([file_to_key(f) => geotiff(SimpleSDMPredictor, joinpath(raster_path, f)) for f in rasters])

# Richness
richness = geotiff(SimpleSDMResponse, joinpath(@__DIR__, "mask.tif"))
richness.grid[findall(!isnothing, richness.grid)] .= 0.0
for (k,v) in ranges
    richness.grid[findall(!isnothing, v.grid)] .+= 1.0
end