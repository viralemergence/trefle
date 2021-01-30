clover_df = "clover.csv" |> CSV.File |> DataFrame
hosts = sort(unique(clover_df.Host))
viruses = sort(unique(clover_df.Virus))
A = zeros(Bool, (length(viruses), length(hosts)))
clover = BipartiteNetwork(A, viruses, hosts)
for clover_row in eachrow(clover_df)
    clover[clover_row.Virus, clover_row.Host] = true
end