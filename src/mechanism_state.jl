# State information pertaining to a single joint
type JointState{X<:Real, M<:Real, C<:Real}
    joint::Joint{M}
    q::VectorSegment{X}
    v::VectorSegment{X}
    beforeJointToParent::Transform3D{C}
    afterJointToParent::Transform3D{C}
    twist::Twist{C}
    biasAcceleration::SpatialAcceleration{C}
    motionSubspace:: JointGeometricJacobian{C}

    function JointState(joint::Joint{M}, beforeJointToParent::Transform3D{C}, q::VectorSegment{X}, v::VectorSegment{X})
        new(joint, q, v, beforeJointToParent)
    end
end

function JointState{X, M}(joint::Joint{M}, beforeJointToParent::Transform3D{M}, q::VectorSegment{X}, v::VectorSegment{X})
    C = promote_type(M, X)
    JointState{X, M, C}(joint, convert(Transform3D{C}, beforeJointToParent), q, v)
end

configuration(state::JointState) = state.q
velocity(state::JointState) = state.v
configuration_range(state::JointState) = first(parentindexes(state.q))
velocity_range(state::JointState) = first(parentindexes(state.v))
parent_frame(state::JointState) = state.beforeJointToParent.to
transform(state::JointState) = state.beforeJointToParent * joint_transform(state.joint, state.q) # FIXME
twist(state::JointState) = change_base(joint_twist(state.joint, state.q, state.v), parent_frame(state))
bias_acceleration(state::JointState) = change_base(bias_acceleration(state.joint, state.q, state.v), parent_frame(state))
motion_subspace(state::JointState) = change_base(motion_subspace(state.joint, state.q), parent_frame(state))
zero_configuration!(state::JointState) = (zero_configuration!(state.joint, state.q))
rand_configuration!(state::JointState) = (rand_configuration!(state.joint, state.q))

# State information pertaining to a single rigid body
type RigidBodyState{M<:Real, C<:Real}
    body::RigidBody{M}
    transformToWorld::Transform3D{C}
    twist::Twist{C}
    biasAcceleration::SpatialAcceleration{C}
    motionSubspace:: JointGeometricJacobian{C} # in world frame
    inertia::SpatialInertia{C}
    crbInertia::SpatialInertia{C}

    function RigidBodyState(body::RigidBody{M}, isroot::Bool)
        ret = new(body)
        if isroot
            ret.transformToWorld = Transform3D(C, body.frame)
            ret.twist = zero(Twist{C}, body.frame, body.frame, body.frame)
            ret.biasAcceleration = zero(SpatialAcceleration{C}, body.frame, body.frame, body.frame)
        end
        ret
    end
end

RigidBodyState{M, X}(body::RigidBody{M}, ::Type{X}, isroot::Bool) = RigidBodyState{M, promote_type(M, X)}(body, isroot)

# MechanismStateCacheStatus stores information regarding what cached variables in a MechanismState are up to date.
immutable MechanismStateCacheStatus
    transforms::Bool
    twists::Bool
    biasAccelerations::Bool
    motionSubspaces::Bool
    inertias::Bool
    crbInertias::Bool

    MechanismStateCacheStatus() = new(false, false, false, false, false, false)

    function MechanismStateCacheStatus(transforms::Bool, twists::Bool, biasAccelerations::Bool, motionSubspaces::Bool, inertias::Bool, crbInertias::Bool)
        inertias = inertias || crbInertias
        twists = twists || biasAccelerations
        transforms = transforms || twists || biasAccelerations || motionSubspaces || inertias
        new(transforms, twists, biasAccelerations, motionSubspaces, inertias, crbInertias)
    end
end

# Note: this constructor allocates and is only meant to be used during codegen.
function MechanismStateCacheStatus{N}(request::NTuple{N, Symbol})
    allowedValues = fieldnames(MechanismStateCacheStatus)
    for element ∈ request
        if element ∉ allowedValues
            error("Unrecognized cache status field: $(string(element)).")
        end
    end

    # TODO: macro
    transforms = :transforms ∈ request
    twists = :twists ∈ request
    biasAccelerations = :biasAccelerations ∈ request
    motionSubspaces = :motionSubspaces ∈ request
    inertias = :inertias ∈ request
    crbInertias = :crbInertias ∈ request
    MechanismStateCacheStatus(transforms, twists, biasAccelerations, motionSubspaces, inertias, crbInertias)
end

# State of an entire Mechanism
# The state information pertaining to rigid bodies in toposortedStateVertices is currently stored in the Mechanism's root frame.
type MechanismState{X<:Real, M<:Real, C<:Real}
    mechanism::Mechanism{M}
    q::Vector{X}
    v::Vector{X}
    toposortedStateVertices::Vector{TreeVertex{RigidBodyState{M, C}, JointState{X, M, C}}}
    nonRootStateVertices::VectorSegment{TreeVertex{RigidBodyState{M, C}, JointState{X, M, C}}}
    cacheStatus::MechanismStateCacheStatus

    function MechanismState(::Type{X}, mechanism::Mechanism{M})
        q = Vector{X}(num_positions(mechanism))
        v = zeros(X, num_velocities(mechanism))
        rootBodyState = RigidBodyState(root_body(mechanism), X, true)
        tree = Tree{RigidBodyState{M, C}, JointState{X, M, C}}(rootBodyState)
        jointStates = Dict{Joint{M}, JointState{X, M, C}}()
        bodyStates = Dict{RigidBody{M}, RigidBodyState{M, C}}()
        for vertex in filter(x -> !isroot(x), mechanism.toposortedTree)
            body = vertex_data(vertex)
            bodyState = RigidBodyState(body, X, false)
            joint = edge_to_parent_data(vertex)
            parentBody = vertex_data(parent(vertex))
            parentStateVertex = findfirst(v -> vertex_data(v).body == parentBody, tree)
            qJoint = view(q, mechanism.qRanges[joint])
            vJoint = view(v, mechanism.vRanges[joint])
            beforeJointToParent = mechanism.jointToJointTransforms[joint]
            jointState = JointState(joint, beforeJointToParent, qJoint, vJoint)
            insert!(parentStateVertex, bodyState, jointState)
            zero_configuration!(joint, qJoint)
        end
        vertices = toposort(tree)
        new(mechanism, q, v, vertices, view(vertices, 2 : length(vertices)), MechanismStateCacheStatus())
    end
end
MechanismState{X, M}(t::Type{X}, mechanism::Mechanism{M}) = MechanismState{X, M, promote_type(X, M)}(t, mechanism)

show{X, M, C}(io::IO, ::MechanismState{X, M, C}) = print(io, "MechanismState{$X, $M, $C}(…)")
num_positions(state::MechanismState) = length(state.q)
num_velocities(state::MechanismState) = length(state.v)
state_vector_eltype{X, M, C}(state::MechanismState{X, M, C}) = X
mechanism_eltype{X, M, C}(state::MechanismState{X, M, C}) = M
cache_eltype{X, M, C}(state::MechanismState{X, M, C}) = C
state_vertex(state::MechanismState, body::RigidBody) = state.toposortedStateVertices[findfirst(v -> vertex_data(v).body == body, state.toposortedStateVertices)] # FIXME: linear time
state_vertex(state::MechanismState, joint::Joint) = state.toposortedStateVertices[findfirst(v -> !isroot(v) && edge_to_parent_data(v).joint == joint, state.toposortedStateVertices)] # FIXME: linear time
configuration(state::MechanismState, joint::Joint) = edge_to_parent_data(state_vertex(state, joint)).q
velocity(state::MechanismState, joint::Joint) = edge_to_parent_data(state_vertex(state, joint)).v
non_root_vertices(state::MechanismState) = state.nonRootStateVertices

function setdirty!(state::MechanismState)
    state.cacheStatus = MechanismStateCacheStatus()
end

function zero_configuration!(state::MechanismState)
    for vertex in non_root_vertices(state)
        zero_configuration!(edge_to_parent_data(vertex))
    end
    setdirty!(state)
end

function zero_velocity!(state::MechanismState)
    X = eltype(state.v)
    fill!(state.v,  zero(X))
    setdirty!(state)
end

zero!(state::MechanismState) = begin zero_configuration!(state); zero_velocity!(state) end

function rand_configuration!(state::MechanismState)
    for vertex in non_root_vertices(state)
        rand_configuration!(edge_to_parent_data(vertex))
    end
    setdirty!(state)
end

function rand_velocity!(state::MechanismState)
    rand!(state.v)
    setdirty!(state)
end

rand!(state::MechanismState) = begin rand_configuration!(state); rand_velocity!(state) end

configuration_vector(state::MechanismState) = state.q
velocity_vector(state::MechanismState) = state.v
state_vector(state::MechanismState) = [configuration_vector(state); velocity_vector(state)]

configuration_vector{T}(state::MechanismState, path::Path{RigidBody{T}, Joint{T}}) = vcat([configuration(state, joint) for joint in path.edgeData]...) # TODO: packing version
velocity_vector{T}(state::MechanismState, path::Path{RigidBody{T}, Joint{T}}) = vcat([velocity(state, joint) for joint in path.edgeData]...) # TODO: packing version

function set_configuration!(state::MechanismState, joint::Joint, q::AbstractVector)
    configuration(state, joint)[:] = q
    setdirty!(state)
end

function set_velocity!(state::MechanismState, joint::Joint, v::AbstractVector)
    velocity(state, joint)[:] = v
    setdirty!(state)
end

function set_configuration!(state::MechanismState, q::AbstractVector)
    copy!(state.q, q)
    setdirty!(state)
end

function set_velocity!(state::MechanismState, v::AbstractVector)
    copy!(state.v, v)
    setdirty!(state)
end

function set!(state::MechanismState, x::AbstractVector)
    nq = num_positions(state)
    nv = num_velocities(state)
    length(x) == nq + nv || error("wrong size")
    @inbounds copy!(state.q, 1, x, 1, nq)
    @inbounds copy!(state.v, 1, x, nq + 1, nv)
    setdirty!(state)
end

@generated function update_cache!{request}(state::MechanismState, ::Type{Val{request}})
    desiredStatus = MechanismStateCacheStatus(request)

    vertexUpdateForwardPass = quote
        parentVertex = parent(vertex)
        jointState = edge_to_parent_data(vertex)
        bodyState = vertex_data(vertex)
    end

    if desiredStatus.transforms
        vertexUpdateForwardPass = quote
            $vertexUpdateForwardPass
            jointState.afterJointToParent = jointState.beforeJointToParent * joint_transform(jointState.joint, jointState.q)
            bodyState.transformToWorld = transform_to_root(parentVertex) * jointState.afterJointToParent
        end
    end

    if desiredStatus.twists
        vertexUpdateForwardPass = quote
            $vertexUpdateForwardPass
            jointState.twist = change_base(joint_twist(jointState.joint, jointState.q, jointState.v), parent_frame(jointState))
            bodyState.twist = twist_wrt_world(parentVertex) + transform(jointState.twist, transform_to_root(vertex))
        end
    end

    if desiredStatus.biasAccelerations
        vertexUpdateForwardPass = quote
            $vertexUpdateForwardPass
            jointState.biasAcceleration = change_base(bias_acceleration(jointState.joint, jointState.q, jointState.v), parent_frame(jointState))
            parentBias = bias_acceleration(parentVertex)
            toRoot = transform_to_root(vertex)
            jointBias = jointState.biasAcceleration
            twistWrtWorld = transform(twist_wrt_world(vertex), inv(toRoot)) # TODO
            jointTwist = twist(jointState)
            jointBias = transform(jointBias, toRoot, twistWrtWorld, jointTwist)
            bodyState.biasAcceleration = parentBias + jointBias
        end
    end

    if desiredStatus.motionSubspaces
        vertexUpdateForwardPass = quote
            $vertexUpdateForwardPass
            jointState.motionSubspace = change_base(motion_subspace(jointState.joint, jointState.q), parent_frame(jointState))
            bodyState.motionSubspace = transform(jointState.motionSubspace, transform_to_root(vertex))
        end
    end

    if desiredStatus.inertias
        vertexUpdateForwardPass = quote
            $vertexUpdateForwardPass
            bodyState.inertia = transform(spatial_inertia(bodyState.body), transform_to_root(vertex))
        end
    end

    if desiredStatus.crbInertias
        vertexUpdateForwardPass = quote
            # initialize crb inertias; finish computation in reverse pass
            $vertexUpdateForwardPass
            bodyState.crbInertia = bodyState.inertia
        end
    end

    ret = quote
        state.cacheStatus = $desiredStatus # do this first so subsequent that cache checks pass
        nonRootVertices = non_root_vertices(state)
        for vertex in nonRootVertices
            $vertexUpdateForwardPass
        end
    end

    if desiredStatus.crbInertias
        ret = quote
            $ret
            # reverse pass for crb inertias
            for i = length(nonRootVertices) : -1 : 1
                vertex = nonRootVertices[i]
                vertex_data(parent(vertex)).crbInertia += crb_inertia(vertex)
            end
        end
    end
    return ret
end

# the following functions return quantities expressed in world frame and w.r.t. world frame (where applicable)
# they do not perform any cache status checks!
@inline function transform_to_root{X, M, C}(vertex::TreeVertex{RigidBodyState{M, C}, JointState{X, M, C}})
    vertex_data(vertex).transformToWorld
end

@inline function twist_wrt_world{X, M, C}(vertex::TreeVertex{RigidBodyState{M, C}, JointState{X, M, C}})
    vertex_data(vertex).twist
end

@inline function bias_acceleration{X, M, C}(vertex::TreeVertex{RigidBodyState{M, C}, JointState{X, M, C}})
    vertex_data(vertex).biasAcceleration
end

@inline function motion_subspace{X, M, C}(vertex::TreeVertex{RigidBodyState{M, C}, JointState{X, M, C}})
    vertex_data(vertex).motionSubspace
end

@inline function spatial_inertia{X, M, C}(vertex::TreeVertex{RigidBodyState{M, C}, JointState{X, M, C}})
    vertex_data(vertex).inertia
end

@inline function crb_inertia{X, M, C}(vertex::TreeVertex{RigidBodyState{M, C}, JointState{X, M, C}})
    vertex_data(vertex).crbInertia
end

function newton_euler{X, M, C}(vertex::TreeVertex{RigidBodyState{M, C}, JointState{X, M, C}}, accel::SpatialAcceleration)
    inertia = spatial_inertia(vertex)
    twist = twist_wrt_world(vertex)
    newton_euler(inertia, accel, twist)
end

momentum{X, M, C}(vertex::TreeVertex{RigidBodyState{M, C}, JointState{X, M, C}}) = spatial_inertia(vertex) * twist_wrt_world(vertex)
momentum_rate_bias{X, M, C}(vertex::TreeVertex{RigidBodyState{M, C}, JointState{X, M, C}}) = newton_euler(vertex, bias_acceleration(vertex))
kinetic_energy{X, M, C}(vertex::TreeVertex{RigidBodyState{M, C}, JointState{X, M, C}}) = kinetic_energy(spatial_inertia(vertex), twist_wrt_world(vertex))

function configuration_derivative!{X}(out::AbstractVector{X}, state::MechanismState{X})
    for vertex in non_root_vertices(state)
        jointState = edge_to_parent_data(vertex)
        q = configuration(jointState)
        v = velocity(jointState)
        q̇ = UnsafeVectorView(out, configuration_range(jointState))
        velocity_to_configuration_derivative!(jointState.joint, q̇, q, v)
    end
end

function configuration_derivative{X}(state::MechanismState{X})
    ret = Vector{X}(num_positions(state.mechanism))
    configuration_derivative!(ret, state)
    ret
end

function transform_to_root(state::MechanismState, frame::CartesianFrame3D)
    body = state.mechanism.bodyFixedFrameToBody[frame]
    tf = transform_to_root(state_vertex(state, body))
    if tf.from != frame
        tf = tf * find_body_fixed_frame_definition(state.mechanism, body, frame) # TODO: consider caching
    end
    tf
end

motion_subspace(state::MechanismState, joint::Joint) = motion_subspace(state_vertex(state, joint))

for fun in (:twist_wrt_world, :bias_acceleration, :spatial_inertia, :crb_inertia, :momentum, :momentum_rate_bias, :kinetic_energy)
    @eval $fun{X, M, C}(state::MechanismState{X, M, C}, body::RigidBody{M}) = $fun(state_vertex(state, body))
end

for fun in (:momentum, :momentum_rate_bias, :kinetic_energy)
    @eval $fun{X, M, C}(state::MechanismState{X, M, C}) = sum($fun, non_root_vertices(state))
    @eval $fun{X, M, C}(state::MechanismState{X, M, C}, body_itr) = sum($fun(state, body) for body in body_itr)
end

function relative_transform(state::MechanismState, from::CartesianFrame3D, to::CartesianFrame3D)
    inv(transform_to_root(state, to)) * transform_to_root(state, from)
end

function relative_twist(state::MechanismState, body::RigidBody, base::RigidBody)
    -twist_wrt_world(state, base) + twist_wrt_world(state, body)
 end

function relative_twist(state::MechanismState, bodyFrame::CartesianFrame3D, baseFrame::CartesianFrame3D)
    twist = relative_twist(state, state.mechanism.bodyFixedFrameToBody[bodyFrame], state.mechanism.bodyFixedFrameToBody[baseFrame])
    Twist(bodyFrame, baseFrame, twist.frame, twist.angular, twist.linear)
end

for VectorType in (:Point3D, :FreeVector3D, :Twist, :Momentum, :Wrench)
    @eval begin
        function transform(state::MechanismState, v::$VectorType, to::CartesianFrame3D)::similar_type(typeof(v), promote_type(cache_eltype(state), eltype(v)))
            # TODO: consider transforming in steps, so that computing the relative transform is not necessary
            v.frame == to ? v : transform(v, relative_transform(state, v.frame, to))
        end
    end
end

function transform(state::MechanismState, accel::SpatialAcceleration, to::CartesianFrame3D)
    accel.frame == to && return accel # nothing to be done
    oldToRoot = transform_to_root(state, accel.frame)
    rootToOld = inv(oldToRoot)
    twistOfBodyWrtBase = transform(relative_twist(state, accel.body, accel.base), rootToOld)
    twistOfOldWrtNew = transform(relative_twist(state, accel.frame, to), rootToOld)
    oldToNew = inv(transform_to_root(state, to)) * oldToRoot
    transform(accel, oldToNew, twistOfOldWrtNew, twistOfBodyWrtBase)
end

function local_coordinates!(state::MechanismState, ϕ::StridedVector, ϕd::StridedVector, q0::StridedVector)
    mechanism = state.mechanism
    for vertex in non_root_vertices(state)
        jointState = edge_to_parent_data(vertex)
        qRange = configuration_range(jointState)
        vRange = velocity_range(jointState)
        ϕjoint = UnsafeVectorView(ϕ, vRange)
        ϕdjoint = UnsafeVectorView(ϕd, vRange)
        q0joint = UnsafeVectorView(q0, qRange)
        qjoint = configuration(jointState)
        vjoint = velocity(jointState)
        local_coordinates!(jointState.joint, ϕjoint, ϕdjoint, q0joint, qjoint, vjoint)
    end
end

function global_coordinates!(state::MechanismState, q0::StridedVector, ϕ::StridedVector)
    mechanism = state.mechanism
    for vertex in non_root_vertices(state)
        jointState = edge_to_parent_data(vertex)
        q0joint = UnsafeVectorView(q0, configuration_range(jointState))
        ϕjoint = UnsafeVectorView(ϕ, velocity_range(jointState))
        qjoint = configuration(jointState)
        global_coordinates!(jointState.joint, qjoint, q0joint, ϕjoint)
    end
end
