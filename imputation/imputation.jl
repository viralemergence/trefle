using DataFrames
using DataFramesMeta
import CSV
using StatsPlots
using Statistics
using EcologicalNetworks
using EcologicalNetworksPlots
using TSne
using GLM

# Get the main predictions
predictions = DataFrame(CSV.File("hpc/outputs/predictions.csv"))

logistic = (x) -> 1.0 / (1.0 + exp(-x))
logit = (p) -> log(p/(1.0-p))

predictions.evidence = (predictions.updated ./ predictions.initial) .- 1.0
predictions.P = logistic.(predictions.evidence)

# Get clover
clover_df = joinpath("data", "clover.csv") |> CSV.File |> DataFrame
hosts = sort(unique(clover_df.Host))
viruses = sort(unique(clover_df.Virus))
A = zeros(Bool, (length(viruses), length(hosts)))
clover = BipartiteNetwork(A, viruses, hosts)
for clover_row in eachrow(clover_df)
    clover[clover_row.Virus, clover_row.Host] = true
end

# Correct the dataframe
for r in eachrow(predictions)
    r.value = clover[r.virus_pos, r.host_pos]
end

imputed = @linq predictions |>
    where(:value .== false) |>
    where(:P .>= 0.846847) |>
    select(:value, :host, :virus, :evidence, :P) |>
    orderby(:evidence)

CSV.write("artifacts/imputed_associations.csv", select(imputed, Not(:value)))

zoonoses = @where(imputed, :host .== "Homo sapiens")
CSV.write("artifacts/zoonoses.csv", select(zoonoses, Not(:value)))

imputed_clover = copy(clover)
for r in eachrow(imputed)
    imputed_clover[r.virus, r.host] = true
end

# Embeddings for Nardus
nohuman = filter(!isequal("Homo sapiens"), species(clover, dims=2))
Lc, Rc = rdpg(clover, 12)
Lt, Rt = rdpg(imputed_clover, 12)
Lnc, Rnc = rdpg(clover[:,nohuman], 12)
Lnt, Rnt = rdpg(imputed_clover[:,nohuman], 12)
df = DataFrame(virus = String[], rank = Int[], clover = Float64[], trefle = Float64[], clover_nohuman = Float64[], trefle_nohuman = Float64[])
for i in 1:size(Lc, 1)
    for j in 1:size(Lc, 2)
        push!(df, (species(clover; dims=1)[i], j, Lc[i,j], Lt[i,j], Lnc[i,j], Lnt[i,j]))
    end
end
CSV.write(joinpath("artifacts", "viral_subspace.csv"), df)


# Writing the predictions
d = DataFrame(virus = String[], host = String[])
for int in interactions(imputed_clover)
    push!(d, (int.from, int.to))
end
sort!(d, :virus)
CSV.write("artifacts/trefle.csv", d)

# Degree distribution
function Pk(N::T; dims::Union{Nothing,Integer}=nothing) where {T<:AbstractEcologicalNetwork}
    deg = collect(values(degree(N; dims=dims)))
    u = sort(unique(deg))
    p = zeros(Float64, length(u))
    for (i,x) in enumerate(u)
        p[i] = sum(deg.==x)
    end
    return (u, p./sum(p))
end

#=
plot(Pk(clover; dims=1)..., lab="Raw data", c=:grey, ls=:dash)
plot!(Pk(imputed_clover; dims=1)..., lab="Post imputation", c=:black)
xaxis!(:log, "Virus degree", (1, 1e3))
yaxis!(:log, (0.0005, 1.0), "Frequency")
savefig("mainfigs/degree-virus.png")

plot(Pk(clover; dims=2)..., lab="Raw data", c=:grey, ls=:dash)
plot!(Pk(imputed_clover; dims=2)..., lab="Post imputation", c=:black)
xaxis!(:log, "Host degree", (1, 1e3))
yaxis!(:log, (0.0005, 1.0), "Frequency")
savefig("mainfigs/degree-host.png")

plot(Pk(clover)..., lab="Raw data", c=:grey, ls=:dash)
plot!(Pk(imputed_clover)..., lab="Post imputation", c=:black)
xaxis!(:log, "Degree", (1, 1e3))
yaxis!(:log, (0.0005, 1.0), "Frequency")
savefig("mainfigs/degree-global.png")
=#

# Embedding
UCLOV = EcologicalNetworks.mirror(convert(UnipartiteNetwork, clover))
emb_clover = tsne(convert.(Float64, Array(UCLOV.edges)), 2, 0, 2000, 6)

IO = initial(RandomInitialLayout, UCLOV)
for (i,s) in enumerate(species(UCLOV))
    IO[s].x = emb_clover[i,1]
    IO[s].y = emb_clover[i,2]
end

UIMPT = EcologicalNetworks.mirror(convert(UnipartiteNetwork, imputed_clover))
emb_impt = tsne(convert.(Float64, Array(UIMPT.edges)), 2, 0, 2000, 6)

IM = initial(RandomInitialLayout, UIMPT)
for (i,s) in enumerate(species(UIMPT))
    IM[s].x = emb_impt[i,1]
    IM[s].y = emb_impt[i,2]
end

p_orig = scatter(IO, clover, bipartite=true, nodesize=degree(clover), msc=:skyblue, msw=0.5, aspectratio=1, dpi=300)
p_imput = scatter(IM, imputed_clover, bipartite=true, nodesize=degree(imputed_clover), msc=:teal, msw=0.5, aspectratio=1, dpi=300)

plot(p_orig, p_imput)
savefig("figures/before-after.png")

# Number of paths
top_zoo = zoonoses.virus#vec(zoonoses[end-19:end,:].virus)
TMP = UnipartiteQuantitativeNetwork(number_of_paths(UCLOV; n=3), EcologicalNetworks._species_objects(UCLOV)...)
i = interactions(TMP)
filter!(int -> int.from != int.to, i)
filter!(int -> "Homo sapiens" in [int.from], i)
filter!(int -> int.to in top_zoo, i)
df = DataFrame(virus = String[], paths = Int64[])
for int in i
    push!(df, (int.to, int.strength))
end

path_zoo = leftjoin(zoonoses, df, on=:virus)

# Residuals from number of paths
X = log.(path_zoo.evidence)
y = log.(path_zoo.paths)

regmodel = lm(permutedims(permutedims(X)), y)

path_zoo.residuals = residuals(regmodel)
sort!(path_zoo, :residuals, rev=true)
select!(path_zoo, Not(:value))
select!(path_zoo, Not(:host))


@df path_zoo scatter(:evidence, :paths, frame=:box, lab="", dpi=400, marker_z=path_zoo.residuals, aspectratio=1, clim=(-3,3), c=:PRGn)
xaxis!(:log, (0.9, 500), "Evidence of interaction")
yaxis!(:log, (0.9, 500), "Number of paths to human")
plot!(x -> coef(regmodel)[1]x, lab="", lw=2.0, c=:darkgrey, ls=:dash)
savefig("figures/number_of_paths.png")

tax = unique(select(clover_df, [:VirusClass, :VirusOrder, :VirusFamily, :VirusGenus, :Virus]))
path_zoo = leftjoin(path_zoo, tax, on = :virus => :Virus, matchmissing = :equal)

@df dropmissing(select(path_zoo, [:residuals, :VirusOrder])) dotplot(:VirusOrder, :residuals, group=:VirusOrder, legend=:outerleft, size=(1000, 400), frame=:zerolines)
xaxis!(rotation=90)

CSV.write(joinpath("artifacts", "number_of_paths.csv"), path_zoo)


# Overlap analysis
overlap_results = DataFrame(sp = String[], step = Symbol[], score = Float64[])

raw_ajs = AJS(clover; dims=2)
raw_human = filter(p -> "Homo sapiens" in p.first, raw_ajs)
for r in raw_human
    sp = filter(p -> p != "Homo sapiens", collect(r.first))[1]
    push!(overlap_results, (sp, :initial, r.second))
end

imp_ajs = AJS(imputed_clover; dims=2)
imp_human = filter(p -> "Homo sapiens" in p.first, imp_ajs)
for r in imp_human
    sp = filter(p -> p != "Homo sapiens", collect(r.first))[1]
    push!(overlap_results, (sp, :imputed, r.second))
end

top_10_initial = @linq overlap_results |>
    where(:step .== ^(:initial)) |>
    orderby(:score)
top_10_initial = last(top_10_initial, 10)

top_10_imputed = @linq overlap_results |>
    where(:step .== ^(:imputed)) |>
    orderby(:score)
top_10_imputed = last(top_10_imputed, 10)

plot(grid=false, xticks=((1,2),["Initial", "Imputed"]), legend=false, dpi=300)
raw_order = LinRange(minimum(top_10_initial.score)-0.02, maximum(top_10_initial.score)+0.02, 10)
for (i,r) in enumerate(raw_order)
    plot!([0.88, 1.0], [raw_order[i], top_10_initial.score[i]], lab="", c=:darkgrey)
end
annotate!([(0.85, raw_order[i], text(r.sp, :right, 6)) for (i,r) in enumerate(eachrow(top_10_initial))])
imp_order = LinRange(minimum(top_10_imputed.score)-0.02, maximum(top_10_imputed.score)+0.02, 10)
for (i,r) in enumerate(imp_order)
    plot!([2.0, 2.12], [top_10_imputed.score[i], r], lab="", c=:darkgrey)
end
annotate!([(2.15, imp_order[i], text(r.sp, :left, 6)) for (i,r) in enumerate(eachrow(top_10_imputed))])
for rinit in eachrow(top_10_initial)
    for rimpt in eachrow(top_10_imputed)
        if rinit.sp == rimpt.sp
            plot!([1, 2], [rinit.score, rimpt.score], c=:grey, ls=:dot)
        end
    end
end
@df top_10_initial scatter!(fill(1, 10), :score, c=:black, msw=0.0)
@df top_10_imputed scatter!(fill(2, 10), :score, c=:black, msw=0.0)
xaxis!((-0.1,3.1), false)
yaxis!((0.0,0.75), "Similarity")
savefig("figures/human-similarity.png")


# Co-occurrences in datasets
imputed.cooc = fill(false, size(imputed, 1))
for i in 1:size(imputed, 1)
    v = imputed.virus[i]
    h = imputed.host[i]
    dbhost = unique(@where(clover_df, :Host .== h).Database)
    dbvirus = unique(@where(clover_df, :Virus .== v).Database)
    imputed.cooc[i] = length(dbvirus ∩ dbhost) > 0
end

@df imputed density(:P, group=:cooc, frame=:box, legend=:topleft, dpi=300)
xaxis!((0.846847, 1.0), "Probability")
yaxis!((0, 45),)
savefig("figures/probability-by-cooccurrence.png")


# Phylogeny - rank correlation
phylodist = DataFrame(CSV.File("artifacts/phylo_distance_to_human.csv"))
rename!(phylodist, "Column1" => "sp")
phylodist.sp = map(n -> replace(n, "_" => " "), phylodist.sp)

sharing_results = DataFrame(sp = String[], step = Symbol[], count = Int64[])

raw_share = overlap(clover; dims=2)
for r in filter(p -> "Homo sapiens" in p.first, raw_share)
    sp = filter(p -> p != "Homo sapiens", collect(r.first))[1]
    push!(sharing_results, (sp, :initial, r.second))
end

imp_share = overlap(imputed_clover; dims=2)
for r in filter(p -> "Homo sapiens" in p.first, imp_share)
    sp = filter(p -> p != "Homo sapiens", collect(r.first))[1]
    push!(sharing_results, (sp, :imputed, r.second))
end

sharing_results = leftjoin(sharing_results, select(overlap_results, :sp, :score), on=:sp)
phyloverlap = leftjoin(sharing_results, phylodist, on=:sp)
CSV.write("artifacts/sharing-phylogeny.csv", phyloverlap)
