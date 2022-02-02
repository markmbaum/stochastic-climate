using DrWatson
@quickactivate "Stochastic Climate"
push!(LOAD_PATH, srcdir())
using StochasticClimate
using PyPlot
using Base.Threads: nthreads, @threads

pygui(true)

##

#number of unique simulations
N = 200*nthreads()
println("$N ensemble members")

#start time
t₁ = 4.0

#noise magnitude
g = 1.5e3

##

#vector of vectors
U = Vector{Vector{Float64}}(undef, N)
@threads for i ∈ 1:N
    #capture only the temperature results
    U[i] = simulate(initparams(t₁=t₁, g=g))[3]
end
#combine all the resutls
T = vcat(U...)

##

c = 100*round(count(x->!isnan(x) & (x < 280), T)/length(T), sigdigits=3)
println("$c % of time in snowball regime")

fig = figure()
hist(T, density=true, log=true, color="gray")
xlabel("Temperature [K]")
ylabel("Density")
title("$c % of time in snowball regime (<280 K)")
fig[:tight_layout]()
