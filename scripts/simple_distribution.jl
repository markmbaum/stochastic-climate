using DrWatson
@quickactivate "Stochastic Climate"
push!(LOAD_PATH, srcdir())
using StochasticClimate
using PyPlot

pygui(true)

##

function Tdistribution(σ, t=𝐭)
    #a bunch of random numbers, centered on reference fCO2
    fCO2 = filter(x->x > 0, σ*randn(1_000_000) .+ fCO2ᵣ)
    #run fCO2 through the simple temperature forcing
    𝒻T.(t, fCO2)
end

##

fig, axs = subplots(3, 1, figsize=(6,7))
for (i,σ) ∈ enumerate([20, 50, 80])
    axs[i][:hist](
        Tdistribution(σ),
        density=true,
        log=true,
        bins=40,
        color="gray"
    )
    axs[i][:set_title]("σ = $σ")
    if i < length(axs)
        axs[i][:set_xticklabels]([])
    end
end
xl = (
    minimum(ax->ax[:get_xlim]()[1], axs),
    maximum(ax->ax[:get_xlim]()[2], axs)
)
for ax ∈ axs
    ax[:set_xlim](xl)
end
fig[:savefig](plotsdir("T_distributions"), dpi=400)