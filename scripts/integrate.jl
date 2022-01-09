using DrWatson
@quickactivate "Stochastic Climate"
push!(LOAD_PATH, srcdir())
using StochasticClimate
using PyPlot

pygui(true)

##

t, fCO2, T = integrate(
    initparams(
        t₁=3.5, #starting time [Gya]
        g=1e4, #magnitude of noise
        τ=1e-3 #strength of "weathering" feedback (smaller is stronger feedback) [Gyr]
    )
)

println("minimum fCO2 = $(minimum(fCO2)) ppm")

fig = figure()

subplot(3,1,1)
plot(𝐭 .- t, map(x -> x <= 0 ? NaN : log10(x/fCO2₀), fCO2))
ylabel("fCO2 [ppm]")
gca()[:invert_xaxis]()

subplot(3,1,2)
plot(𝐭 .- t, T)
ylabel("Temperature [K]")
xlabel("Time [Gya]")
gca()[:invert_xaxis]()

subplot(3,1,3)
hist(T[(@. !isnan(T))], density=true, log=true, color="gray")
xlabel("Temperature (excluding NaNs) [K]")
ylabel("Density")

fig[:tight_layout]()
