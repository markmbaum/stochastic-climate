module StochasticClimate

using UnPack
using Roots
using BasicInterpolators: ChebyshevInterpolator
using StochasticDiffEq

#------------------------------------------------------------------------------
# physical constants

export F₀, fCO2₀, T₀, A₀, OLR₀, 𝐭

#reference instellation (solar "constant") [W/m^2]
const F₀ = 1366.0
#reference molar concentration of CO2 [ppm]
const fCO2₀ = 280
#reference temperature [K]
const T₀ = 288.0
#reference albedo [-]
const A₀ = 0.3
#reference OLR [W/m^2]
const OLR₀ = (1 - A₀)*F₀/4
#OLR response to temperature
const a = 2.0
#OLR response to pCO2
const b = 5.35
#maximum age (of solar system)
const 𝐭 = 4.5

#------------------------------------------------------------------------------
#component physical equations

export 𝒻α, 𝒻F, 𝒻S, 𝒻OLR, 𝒻A, 𝒻T, 𝒻Tsafe, imbalance

#stellar luminosity fraction over time [Gya]
𝒻α(t) = 1/(1 + (2/5)*(1 - t/𝐭))

#instellation over time [W/m^2]
𝒻F(t) = 𝒻α(t)*F₀

#absorbed radiation [W/m^2]
𝒻S(t, A) = (1 - A)*𝒻F(t)/4

#outgoing longwave radiation [W/m^2]
𝒻OLR(T, fCO2) = OLR₀ + a*(T - T₀) - b*log(fCO2/fCO2₀)

#instantaneous temperature [K]
𝒻T(t, fCO2) = (𝒻S(t, A₀) - OLR₀ + b*log(fCO2/fCO2₀))/a + T₀
𝒻Tsafe(t, fCO2) = fCO2 <= 0 ? NaN : 𝒻T(t, fCO2)

#radiative imbalance [W/m^2]
imbalance(t, T, fCO2) = 𝒻S(t, A₀) - 𝒻OLR(T, fCO2)

#------------------------------------------------------------------------------
# equilibrium fCO2 over time

export 𝒻fCO2ₑ, Χ, tsat

#find the fCO2 value which gives T₀ at some time
𝒻fCO2ₑ(t, Tₑ=T₀) = exp10(find_zero(x -> imbalance(t, Tₑ, exp10(x)), (-12, 12)))

#simple struct to rapidly interpolate χ values instead of root finding
struct Χ{𝒯<:ChebyshevInterpolator} #capital Chi here
    interpolator::𝒯
end

Χ(Tₑ::Real=T₀; N::Int=5) = Χ(ChebyshevInterpolator(t -> log(𝒻fCO2ₑ(t,Float64(Tₑ))), 0.0, 𝐭, N))

(χ::Χ)(t) = exp(χ.interpolator(t))

tsat(χ::Χ) = find_zero(t -> χ(t) - 1e6, (0.0, 𝐭))

#------------------------------------------------------------------------------
# setup & integration

export initparams, integrate

g(u, p, t) = p.g

function 𝒹fCO2(fCO2, p, t)
    #unpack parameters
    @unpack 𝒻χ, τ = p
    #equilibrium 𝒻CO2
    χ = 𝒻χ(t)
    #change in fCO2
    -(fCO2 - χ)/τ
end

function initparams(;
                    Tₑ=288.0, #equilibrium temperature
                    τ=1e-2, #weathering feedback time scale
                    g=1e3, #noise strength
                    dt=1e-1, #time step
                    t₁=2.5, #initial time
                    t₂=4.5, #final time
                    enforcepos=true #whether to include callback preventing negative fCO₂
                    )::NamedTuple
    #named tuple containing integration parameters
    (
        𝒻χ = Χ(Tₑ),
        τ  = Float64(τ),
        g  = Float64(g),
        dt = Float64(dt),
        t₁ = Float64(t₁),
        t₂ = Float64(t₂),
        enforcepos = enforcepos
    )
end

function enforcepositivity()::DiscreteCallback
    DiscreteCallback(
        (u, t, integrator) -> u < 0,
        integrator -> integrator.u = -integrator.u,
        save_positions=(true,true)
    )
end

function integrate(params=initparams())
    @unpack t₁, t₂, 𝒻χ, dt = params
    tspan = (t₁, t₂)
    u₀ = 𝒻χ(t₁)
    #prevent negative fCO₂ or don't
    if params.enforcepos
        println("yes")
        prob = SDEProblem(𝒹fCO2, g, u₀, tspan, params, callback=enforcepositivity())
    else
        prob = SDEProblem(𝒹fCO2, g, u₀, tspan, params)
    end
    sol = solve(prob, SRA3(), dt=dt)
    println(length(sol.t))
    println(minimum(sol.u))
    return sol.t, sol.u, 𝒻Tsafe.(sol.t, sol.u)
end

end