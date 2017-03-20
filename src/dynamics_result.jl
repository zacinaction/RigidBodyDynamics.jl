"""
$(TYPEDEF)

Stores variables related to the dynamics of a `Mechanism`, e.g. the
`Mechanism`'s mass matrix and joint acceleration vector.

Type parameters:
* `M`: the scalar type of the `Mechanism`.
* `T`: the scalar type of the dynamics-related variables.
"""
type DynamicsResult{M<:Number, T<:Number}
    mechanism::Mechanism{M}

    massmatrix::Symmetric{T, Matrix{T}}
    dynamicsbias::Vector{T}
    constraintjacobian::Matrix{T}
    constraintbias::Vector{T}

    v̇::Vector{T}
    λ::Vector{T}
    ṡ::Vector{T}

    # the following are indexed by vertex_index(body). TODO: consider adding a BodyMap type:
    contactwrenches::Vector{Wrench{T}}
    totalwrenches::Vector{Wrench{T}}
    accelerations::Vector{SpatialAcceleration{T}}
    jointwrenches::Vector{Wrench{T}} # TODO: index by joint tree index?
    contact_state_derivatives::Vector{Vector{DefaultSoftContactStateDeriv{T}}}

    # see solve_dynamics! for meaning of the following variables:
    L::Matrix{T} # lower triangular
    A::Matrix{T} # symmetric
    z::Vector{T}
    Y::Matrix{T}

    function (::Type{DynamicsResult{M, T}}){M<:Number, T<:Number}(mechanism::Mechanism{M})
        nq = num_positions(mechanism)
        nv = num_velocities(mechanism)

        nconstraints = mapreduce(num_constraints, +, 0, non_tree_joints(mechanism))

        massmatrix = Symmetric(Matrix{T}(nv, nv), :L)
        dynamicsbias = Vector{T}(nv)
        constraintjacobian = Matrix{T}(nconstraints, nv)
        constraintbias = Vector{T}(nconstraints)

        v̇ = Vector{T}(nv)
        λ = Vector{T}(nconstraints)
        ṡ = Vector{T}(num_additional_states(mechanism))

        contactwrenches = Vector{Wrench{T}}(num_bodies(mechanism))
        totalwrenches = Vector{Wrench{T}}(num_bodies(mechanism))
        accelerations = Vector{SpatialAcceleration{T}}(num_bodies(mechanism))
        jointwrenches = Vector{Wrench{T}}(num_bodies(mechanism))
        contact_state_derivatives = Vector{Vector{DefaultSoftContactStateDeriv{T}}}(num_bodies(mechanism))
        start = 1
        for body in bodies(mechanism)
            derivs = DefaultSoftContactStateDeriv{T}[]
            for point in contact_points(body)
                model = contact_model(point)
                ṡ_part = view(ṡ, start : start + num_states(model) - 1)
                push!(derivs, SoftContactStateDeriv(model, ṡ_part, root_frame(mechanism)))
                start += num_states(model)
            end
            contact_state_derivatives[vertex_index(body)] = derivs
        end

        L = Matrix{T}(nv, nv)
        A = Matrix{T}(nconstraints, nconstraints)
        z = Vector{T}(nv)
        Y = Matrix{T}(nconstraints, nv)

        new{M, T}(mechanism, massmatrix, dynamicsbias, constraintjacobian, constraintbias,
            v̇, λ, ṡ, contactwrenches, totalwrenches, accelerations, jointwrenches, contact_state_derivatives,
            L, A, z, Y)
    end
end

DynamicsResult{M, T}(::Type{T}, mechanism::Mechanism{M}) = DynamicsResult{M, T}(mechanism)

contact_state_derivatives(result::DynamicsResult, body::RigidBody) = result.contact_state_derivatives[vertex_index(body)]
contact_wrench(result::DynamicsResult, body::RigidBody) = result.contactwrenches[vertex_index(body)]
set_contact_wrench!(result::DynamicsResult, body::RigidBody, wrench::Wrench) = (result.contactwrenches[vertex_index(body)] = wrench)
acceleration(result::DynamicsResult, body::RigidBody) = result.accelerations[vertex_index(body)]
set_acceleration!(result::DynamicsResult, body::RigidBody, accel::SpatialAcceleration) = (result.accelerations[vertex_index(body)] = accel)
joint_wrench(result::DynamicsResult, body::RigidBody) = result.jointwrenches[vertex_index(body)]
set_joint_wrench!(result::DynamicsResult, body::RigidBody, wrench::Wrench) = (result.jointwrenches[vertex_index(body)] = wrench)
