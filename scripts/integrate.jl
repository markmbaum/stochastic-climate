using DrWatson
@quickactivate "Stochastic Climate"
push!(LOAD_PATH, srcdir())
using StochasticClimate
using PyPlot
using Base.Threads: nthreads, @threads

pygui(true)

##

params = initparams(
    t₁=2.0, #starting time [Gya]
    g=1e4, #magnitude of noise
    τ=2.5e-3, #strength of "weathering" feedback (smaller is stronger feedback) [Gyr]
    enforcepos=true
)

## show a single integration

t, fCO2, T = integrate(params)

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

N = 100*nthreads()
println("$N ensemble members")
Tₐ = Vector{Vector{Float64}}(undef, N)
@threads for i ∈ 1:N
    Tₐ[i] = integrate(params)[3]
end
Tc = vcat(Tₐ...)
T = filter(x->!isnan(x), Tc)

c = 100*round(count(x->!isnan(x) & (x < 280), Tc)/length(Tc), sigdigits=3)
println("$c % of time in snowball regime")
println("$(size(Tc)[1]-size(T)[1]) NaN Ts in ensemble")

fig = figure()
hist(T, density=true, log=true, color="gray")
xlabel("Temperature (excluding NaNs) [K]")
ylabel("Density")
title("$c % of time in snowball regime (<280 K)")
fig[:tight_layout]()
