using DrWatson
@quickactivate "Stochastic Climate"
push!(LOAD_PATH, srcdir())
using StochasticClimate
using PyPlot
using Base.Threads: nthreads, @threads

pygui(true)

##

params = initparams(
    t₁=4.0, #starting time [Gya]
    g=1e4, #magnitude of noise
    τ=2.5e-3 #strength of "weathering" feedback (smaller is stronger feedback) [Gyr]
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
Tall = Vector{Vector{Float64}}(undef, N)
@threads for i ∈ 1:N
    Tall[i] = integrate(params)[3]
end
Tcat = vcat(Tall...)
T = filter(x->!isnan(x), Tcat)

c = 100*round(count(x->!isnan(x) & (x < 280), Tcat)/length(Tcat), sigdigits=3)
println("$c % of time in snowball regime")

fig = figure()
hist(T, density=true, log=true, color="gray")
xlabel("Temperature (excluding NaNs) [K]")
ylabel("Density")
fig[:tight_layout]()
