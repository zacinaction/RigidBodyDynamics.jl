module Contact

using RigidBodyDynamics
using StaticArrays
using Compat

export NormalForceModel,
    FrictionModel,
    SoftContactModel,
    ContactPoint,
    DefaultContactPoint

export normal_force,
    friction_force,
    reset!,
    dynamics!

export HuntCrossleyModel,
    hunt_crossley_hertz

export
    ViscoelasticCoulombModel,
    ViscoelasticCoulombState

# Types
@compat abstract type NormalForceModel{T} end
normal_force(::NormalForceModel, z, ż) = error("Subtypes must implement")

@compat abstract type FrictionModel{T} end
# TODO: FrictionModel interface

immutable SoftContactModel{N <: NormalForceModel, F <: FrictionModel}
    normal_force_model::N
    friction_model::F
end

type ContactPoint{T, M <: SoftContactModel}
    location::Point3D{SVector{3, T}}
    model::M
end

# Normal contact models
immutable HuntCrossleyModel{T} <: NormalForceModel{T}
    # (2) in Marhefka, Orin, "A Compliant Contact Model with Nonlinear Damping for Simulation of Robotic Systems"
    k::T
    λ::T
    n::T
end

function hunt_crossley_hertz(; k = 50e3, α = 0.2)
    λ = 3/2 * α * k # (12) in Marhefka, Orin
    HuntCrossleyModel(k, λ, 3/2)
end

function normal_force(model::HuntCrossleyModel, z, ż)
    zn = z^model.n
    f = -model.λ * zn * ż - model.k * zn # (2) in Marhefka, Orin
end

# Friction models
immutable ViscoelasticCoulombModel{T} <: FrictionModel{T}
    # See section 11.8 of Featherstone, "Rigid Body Dynamics Algorithms", 2008
    μ::T
    k::T
    b::T
end

# One state for each direction; technically only need one for each tangential direction, but this is easier to work with:
num_states(::ViscoelasticCoulombModel) = 3

type ViscoelasticCoulombState{T, V}
    model::ViscoelasticCoulombModel{T}
    tangential_displacement::FreeVector3D{V}
end

function reset!(state::ViscoelasticCoulombState)
    fill!(state.tangential_displacement.v, 0)
end

function friction_force{T}(state::ViscoelasticCoulombState{T}, fnormal::T, tangential_velocity::FreeVector3D)
    model = state.model
    μ = model.μ
    k = model.k
    b = model.b
    x = convert(FreeVector3D{SVector{3, T}}, state.tangential_displacement)
    v = tangential_velocity

    # compute friction force that would be needed to avoid slip
    fstick = -k * x - b * v

    # limit friction force to lie within friction cone
    fstick_norm² = dot(fstick, fstick)
    fstick_max_norm² = (μ * fnormal)^2
    ftangential = if fstick_norm² > fstick_max_norm²
        fstick * sqrt(fstick_max_norm² / fstick_norm²)
    else
        fstick
    end
end

function dynamics!{T}(ẋ, state::ViscoelasticCoulombState{T}, ftangential::FreeVector3D)
    # TODO: type of ẋ?
    model = state.model
    k = model.k
    b = model.b
    x = state.x
    ẋ .= -(k * x + ftangential) / b
end

@compat const DefaultContactPoint{T} = ContactPoint{T,SoftContactModel{HuntCrossleyModel{T},ViscoelasticCoulombModel{T}}}

end
