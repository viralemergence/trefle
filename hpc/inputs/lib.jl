function lowrank(target, rank)
    f = svd(target)
    f.S[(rank+1):end] .= 0.0
    return f.U * Diagonal(f.S) * f.Vt
end

n1(Y) = fill(mean(Y), size(Y))
n2(Y) = (mean(Y; dims=2).+mean(Y; dims=1))./2.0
n3(Y) = (mean(Y; dims=2).+mean(Y; dims=1).+mean(Y))./3.0

function impute(query, template, position; rank=4)
    initial_value = query[position]
    query[position] = template[position]
    for iteration in 1:20
        _tmp = lowrank(query, rank)
        query[position] = _tmp[position]
    end
    prediction = (position=position, initial=template[position], updated=query[position])
    query[position] = initial_value
    return prediction
end
