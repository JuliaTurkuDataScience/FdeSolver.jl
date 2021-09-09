function α(n, β)

    if n == 0

        return 1 ./ Γ(β .+ 2)

    else

        return ((n .- 1) .^ (β .+ 1) .- 2 .* n .^ (β .+ 1) .+ (n .+ 1) .^ (β .+ 1)) ./ Γ(β .+ 2)

    end

end
