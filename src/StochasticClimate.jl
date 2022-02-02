module StochasticClimate

using UnPack
using Roots
using BasicInterpolators: ChebyshevInterpolator
using StochasticDiffEq

#------------------------------------------------------------------------------
# physical constants

export Fᵣ, fCO2ᵣ, Tᵣ, αᵣ, OLRᵣ, 𝐭

#reference instellation (solar "constant") [W/m^2]
const Fᵣ = 1366.0
#reference molar concentration of CO2 [ppm]
const fCO2ᵣ = 280
#reference temperature [K]
const Tᵣ = 288.0
#reference albedo [-]
const αᵣ = 0.3
#reference OLR [W/m^2]
const OLRᵣ = (1 - αᵣ)*Fᵣ/4
#OLR response to temperature
const a = 2.0
#OLR response to pCO2
const b = 5.35
#maximum age (of solar system)
const 𝐭 = 4.5

#------------------------------------------------------------------------------
#component physical equations

export 𝒻α, 𝒻F, 𝒻S, 𝒻OLR, 𝒻T, 𝒻Tsafe, imbalance

#stellar luminosity fraction over time [Gya]
𝒻α(t=𝐭) = 1/(1 + (2/5)*(1 - t/𝐭))

#instellation over time [W/m^2]
𝒻F(t=𝐭) = 𝒻α(t)*Fᵣ

#absorbed radiation [W/m^2]
𝒻S(t=𝐭, α=αᵣ) = (1 - α)*𝒻F(t)/4

#outgoing longwave radiation [W/m^2]
𝒻OLR(T=Tᵣ, fCO2=fCO2₀) = OLRᵣ + a*(T - Tᵣ) - b*log(fCO2/fCO2ᵣ)

#instantaneous temperature [K]
𝒻T(t=𝐭, fCO2=fCO2₀) = (𝒻S(t, αᵣ) - OLRᵣ + b*log(fCO2/fCO2ᵣ))/a + Tᵣ
𝒻Tsafe(t=𝐭, fCO2=fCO2₀) = fCO2 <= 0 ? NaN : 𝒻T(t, fCO2)

#radiative imbalance [W/m^2]
imbalance(t=𝐭, T=Tᵣ, fCO2=fCO2₀) = 𝒻S(t, αᵣ) - 𝒻OLR(T, fCO2)

#------------------------------------------------------------------------------
# equilibrium fCO2 over time

export 𝒻fCO2ₑ, Χ, tsat

#find the fCO2 value which gives Tᵣ at some time
𝒻fCO2ₑ(t, Tₑ=Tᵣ) = exp10(find_zero(x -> imbalance(t, Tₑ, exp10(x)), (-12, 12)))

#simple struct to rapidly interpolate χ values instead of root finding
struct Χ{𝒯<:ChebyshevInterpolator} #capital Chi here
    interpolator::𝒯
end

Χ(Tₑ::Real=Tᵣ; N::Int=5) = Χ(ChebyshevInterpolator(t -> log(𝒻fCO2ₑ(t,Float64(Tₑ))), 0.0, 𝐭, N))

(χ::Χ)(t) = exp(χ.interpolator(t))

tsat(χ::Χ) = find_zero(t -> χ(t) - 1e6, (0.0, 𝐭))

#------------------------------------------------------------------------------
# setup & integration

export initparams, simulate, integrate

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
                    g=1e2, #noise strength
                    dt=1e-4, #time step
                    t₁=2.5, #initial time
                    t₂=4.5, #final time
                    reflect=true #whether to include callback preventing negative fCO₂
                    )::NamedTuple
    #named tuple containing integration parameters
    (
        𝒻χ = Χ(Tₑ),
        τ  = Float64(τ),
        g  = Float64(g),
        dt = Float64(dt),
        t₁ = Float64(t₁),
        t₂ = Float64(t₂),
        reflect = reflect
    )
end

function reflect!(integrator)::Nothing
    if integrator.u < 0
        integrator.u = -integrator.u
    end
    nothing
end

function reflector()::DiscreteCallback
    # for whatever reason save_positions=(false,true) doesn't
    # accomplish not saving pre-flipped negative CO2
    # so follow Chris Rackauckas' answer here
    # https://stackoverflow.com/questions/69049991/simulating-a-reflecting-boundary-sdeproblem
    # also requires setting save_everystep=false in solve
    DiscreteCallback(
        (u,t,integrator) -> true,
        reflect!,
        save_positions=(false,true)
    )
end

function simulate(params)
    @unpack t₁, t₂, 𝒻χ, dt = params
    tspan = (t₁, t₂)
    u₀ = 𝒻χ(t₁)
    prob = SDEProblem(𝒹fCO2, g, u₀, tspan, params)
    sol = if params.reflect
        solve(
            prob,
            EM(),
            dt=dt,
            callback=reflector(),
            save_everystep=false
        )
    else
        solve(prob, EM(), dt=dt)
    end
    return sol.t, sol.u, 𝒻Tsafe.(sol.t, sol.u)
end

end
