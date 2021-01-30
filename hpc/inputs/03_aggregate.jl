using DataFrames
import CSV
f = readdir()
filter!(n -> startswith(n, "svd-pred"), f)
filter!(n -> endswith(n, "csv"), f)
predict = vcat(DataFrame.(CSV.File.(f))...)
CSV.write("predictions.csv", predict)