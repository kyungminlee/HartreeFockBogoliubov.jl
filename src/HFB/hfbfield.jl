export HFBAmplitude,
       HFBField,
       HFBAmplitudeHint,
       HFBFieldHint

export iscompatible


mutable struct HFBAmplitude
    ρ ::Vector{ComplexF64}
    t ::Vector{ComplexF64}
end

mutable struct HFBField
    Γ ::Vector{ComplexF64}
    Δ ::Vector{ComplexF64}
end

mutable struct HFBAmplitudeHint{O}
    ρ ::Dict{Tuple{O, O, Vector{Int}}, ComplexF64}
    t ::Dict{Tuple{O, O, Vector{Int}}, ComplexF64}
end

mutable struct HFBFieldHint{O}
    Γ ::Dict{Tuple{O, O, Vector{Int}}, ComplexF64}
    Δ ::Dict{Tuple{O, O, Vector{Int}}, ComplexF64}
end


function iscompatible(s1 ::HFBAmplitude, s2::HFBAmplitude) ::Bool
    return (length(s1.ρ) == length(s2.ρ) && length(s1.t) == length(s2.t))
end


function iscompatible(s1 ::HFBField, s2::HFBField) ::Bool
    return (length(s1.Γ) == length(s2.Γ) && length(s1.Δ) == length(s2.Δ))
end


import Base: copy

copy(x::HFBAmplitude) = HFBAmplitude(copy(x.ρ), copy(x.t))
copy(x::HFBField) = HFBField(copy(x.Γ), copy(x.Δ))
copy(x::HFBAmplitudeHint{O}) where {O} = HFBAmplitudeHint{O}(copy(x.ρ), copy(x.t))

import Base.+
import Base.-
import Base.*
import Base./
import Base.\

+(sol1 ::HFBAmplitude) = HFBAmplitude(sol1.ρ, sol1.t)
-(sol1 ::HFBAmplitude) = HFBAmplitude(-sol1.ρ, -sol1.t)
+(sol1 ::HFBAmplitude, sol2 ::HFBAmplitude) = HFBAmplitude(sol1.ρ + sol2.ρ, sol1.t + sol2.t)
-(sol1 ::HFBAmplitude, sol2 ::HFBAmplitude) = HFBAmplitude(sol1.ρ - sol2.ρ, sol1.t - sol2.t)
*(sol1 ::HFBAmplitude, value ::Real) = HFBAmplitude(sol1.ρ * value, sol1.t * value)
*(value ::Real, sol1::HFBAmplitude) = HFBAmplitude(value * sol1.ρ, value * sol1.t)
/(sol1::HFBAmplitude, value ::Real) = HFBAmplitude(sol1.ρ / value, sol1.t / value)
\(value ::Real, sol1::HFBAmplitude) = HFBAmplitude(value \ sol1.ρ, value \ sol1.t)

+(sol1 ::HFBField) = HFBField(sol1.Γ, sol1.Δ)
-(sol1 ::HFBField) = HFBField(-sol1.Γ, -sol1.Δ)
+(sol1 ::HFBField, sol2 ::HFBField) = HFBField(sol1.Γ + sol2.Γ, sol1.Δ + sol2.Δ)
-(sol1 ::HFBField, sol2 ::HFBField) = HFBField(sol1.Γ - sol2.Γ, sol1.Δ - sol2.Δ)
*(sol1 ::HFBField, value ::Real) = HFBField(sol1.Γ * value, sol1.Δ * value)
*(value ::Real, sol1::HFBField) = HFBField(value * sol1.Γ, value * sol1.Δ)
/(sol1::HFBField, value ::Real) = HFBField(sol1.Γ / value, sol1.Δ / value)
\(value ::Real, sol1::HFBField) = HFBField(value \ sol1.Γ, value \ sol1.Δ)

import Base: isapprox

function isapprox(sol1::HFBAmplitude, sol2::HFBAmplitude;
                  atol::Real=sqrt(eps(Float64)),
                  rtol::Real=sqrt(eps(Float64)),
                  nans::Bool=false) ::Bool
    return (isapprox(sol1.ρ, sol2.ρ; rtol=rtol, atol=atol, nans=nans) &&
            isapprox(sol1.t, sol2.t; rtol=rtol, atol=atol, nans=nans))
end

function isapprox(sol1::HFBField, sol2::HFBField;
                  atol::Real=sqrt(eps(Float64)),
                  rtol::Real=sqrt(eps(Float64)),
                  nans::Bool=false) ::Bool
    return (isapprox(sol1.Γ, sol2.Γ; rtol=rtol, atol=atol, nans=nans) &&
            isapprox(sol1.Δ, sol2.Δ; rtol=rtol, atol=atol, nans=nans) )
end
