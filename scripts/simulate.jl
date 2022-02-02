using DrWatson
@quickactivate "Stochastic Climate"
push!(LOAD_PATH, srcdir())
using StochasticClimate
using PyPlot

pygui(true)

##

params = initparams(
    t₁=2.2, #starting time [Gya]
    g=3e3, #magnitude of noise
    τ=2.5e-3, #strength of "weathering" feedback (smaller is stronger feedback) [Gyr]
    reflect=true
)

## show a single integration

t, fCO2, T = simulate(params)

println("minimum fCO2 = $(minimum(fCO2)) ppm")

fig = figure()

subplot(2,1,1)
plot(𝐭 .- t, map(x -> x <= 0 ? NaN : log10(x), fCO2))
ylabel("log₁₀(fCO2) [ppm]")
gca()[:invert_xaxis]()

subplot(2,1,2)
plot(𝐭 .- t, T)
ylabel("Temperature [K]")
xlabel("Time [Gya]")
gca()[:invert_xaxis]()

fig[:tight_layout]()

## temperature statistics for an ensemble of integrations