"""
$(SIGNATURES)

Return the mass of a subtree of a `Mechanism`, rooted at `base` (including
the mass of `base`).
"""
function subtree_mass{T}(base::RigidBody{T}, mechanism::Mechanism{T})
    result = isroot(base, mechanism) ? zero(T) : spatial_inertia(base).mass
    for joint in joints_to_children(base, mechanism)
        result += subtree_mass(successor(joint, mechanism), mechanism)
    end
    result
end

"""
$(SIGNATURES)

Return the total mass of the `Mechanism`.
"""
mass(m::Mechanism) = subtree_mass(root_body(m), m)

"""
$(SIGNATURES)

Compute the center of mass of an iterable subset of a `Mechanism`'s bodies in
the given state.
"""
function center_of_mass{X, M, C}(state::MechanismState{X, M, C}, itr)
    frame = root_frame(state.mechanism)
    com = Point3D(frame, zeros(SVector{3, C}))
    mass = zero(C)
    for body in itr
        inertia = spatial_inertia(body)
        if inertia.mass > zero(C)
            bodyCom = center_of_mass(inertia)
            com += inertia.mass * FreeVector3D(transform(state, bodyCom, frame))
            mass += inertia.mass
        end
    end
    com /= mass
    com
end

"""
$(SIGNATURES)

Compute the center of mass of the whole `Mechanism` in the given state.
"""
center_of_mass(state::MechanismState) = center_of_mass(state, non_root_bodies(state.mechanism))

function _set_jacobian_part!(out::GeometricJacobian, part::GeometricJacobian, nextbaseframe::CartesianFrame3D, startind::Int64)
    @framecheck part.frame out.frame
    @framecheck part.base nextbaseframe
    n = 3 * num_cols(part)
    copy!(out.angular, startind, part.angular, 1, n)
    copy!(out.linear, startind, part.linear, 1, n)
    startind += n
    nextbaseframe = part.body
    startind, nextbaseframe
end

function geometric_jacobian!{A, X, M, C}(out::GeometricJacobian{A}, state::MechanismState{X, M, C},
        path::TreePath{RigidBody{M}, Joint{M}})
    @boundscheck num_velocities(path) == num_cols(out) || error("size mismatch")

    mechanism = state.mechanism
    nextbaseframe = out.base
    startind = 1
    for joint in path.source_to_lca
        S = -motion_subspace_in_world(state, joint)
        startind, nextbaseframe = _set_jacobian_part!(out, S, nextbaseframe, startind)
    end
    for i = length(path.target_to_lca) : -1 : 1
        joint = path.target_to_lca[i]
        S = motion_subspace_in_world(state, joint)
        startind, nextbaseframe = _set_jacobian_part!(out, S, nextbaseframe, startind)
        if i == 1
            @framecheck S.body out.body
        end
    end
    out
end

"""
$(SIGNATURES)

Compute a geometric Jacobian (also known as a basic, or spatial Jacobian)
for a path in the graph of joints and bodies of a `Mechanism`,
in the given state.

See also [`path`](@ref), [`GeometricJacobian`](@ref).
"""
function geometric_jacobian{X, M, C}(state::MechanismState{X, M, C}, path::TreePath{RigidBody{M}, Joint{M}})
    nv = num_velocities(path)
    angular = Matrix{C}(3, nv)
    linear = Matrix{C}(3, nv)
    bodyframe = default_frame(path.target)
    baseframe = default_frame(path.source)
    jac = GeometricJacobian(bodyframe, baseframe, root_frame(state.mechanism), angular, linear)
    geometric_jacobian!(jac, state, path)
end

"""
$(SIGNATURES)

Compute the spatial acceleration of `body` with respect to `base` for
the given state and joint acceleration vector ``\\dot{v}``.
"""
function relative_acceleration(state::MechanismState, body::RigidBody, base::RigidBody, v̇::AbstractVector)
    # This is kind of a strange algorithm due to the availability of bias accelerations w.r.t. world
    # in MechanismState, while computation of the v̇-dependent terms follows the shortest path in the tree.
    # TODO: consider doing everything in body frame

    C = cache_eltype(state)
    mechanism = state.mechanism

    bodyframe = default_frame(body)
    baseframe = default_frame(base)
    rootframe = root_frame(mechanism)

    bodyaccel = zero(SpatialAcceleration{C}, bodyframe, bodyframe, rootframe)
    baseaccel = zero(SpatialAcceleration{C}, baseframe, baseframe, rootframe)

    bias = -bias_acceleration(state, base) + bias_acceleration(state, body)
    while body != base
        do_body = tree_index(body, mechanism) > tree_index(base, mechanism)

        joint = do_body ? joint_to_parent(body, mechanism) : joint_to_parent(base, mechanism)
        S = motion_subspace_in_world(state, joint)
        v̇joint = UnsafeVectorView(v̇, velocity_range(state, joint))
        jointaccel = SpatialAcceleration(S, v̇joint)

        if do_body
            bodyaccel = jointaccel + bodyaccel
            body = predecessor(joint, mechanism)
        else
            baseaccel = jointaccel + baseaccel
            base = predecessor(joint, mechanism)
        end
    end
    (-baseaccel + bodyaccel) + bias
end

"""
$(SIGNATURES)

Compute the gravitational potential energy in the given state, computed as the
negation of the dot product of the gravitational force and the center
of mass expressed in the `Mechanism`'s root frame.
"""
function gravitational_potential_energy{X, M, C}(state::MechanismState{X, M, C})
    m = mass(state.mechanism)
    gravitationalForce = m * state.mechanism.gravitationalAcceleration
    centerOfMass = transform(state, center_of_mass(state), gravitationalForce.frame)
    -dot(gravitationalForce, FreeVector3D(centerOfMass))
 end

 function force_space_matrix_transpose_mul_jacobian!(out::Matrix, rowstart::Int64, colstart::Int64,
        mat::Union{MomentumMatrix, WrenchMatrix}, jac::GeometricJacobian, sign::Int64)
     @framecheck(jac.frame, mat.frame)
     n = num_cols(mat)
     m = num_cols(jac)
     @boundscheck (rowstart > 0 && rowstart + n - 1 <= size(out, 1)) || error("size mismatch")
     @boundscheck (colstart > 0 && colstart + m - 1 <= size(out, 2)) || error("size mismatch")

     # more efficient version of
     # view(out, rowstart : rowstart + n - 1, colstart : colstart + n - 1)[:] = mat.angular' * jac.angular + mat.linear' * jac.linear

     @inbounds begin
         for col = 1 : m
             outcol = colstart + col - 1
             for row = 1 : n
                 outrow = rowstart + row - 1
                 out[outrow, outcol] = zero(eltype(out))
                 for i = 1 : 3
                     out[outrow, outcol] += flipsign(mat.angular[i, row] * jac.angular[i, col], sign)
                     out[outrow, outcol] += flipsign(mat.linear[i, row] * jac.linear[i, col], sign)
                 end
             end
         end
     end
 end

const mass_matrix_doc = """Compute the joint-space mass matrix
(also known as the inertia matrix) of the `Mechanism` in the given state, i.e.,
the matrix ``M(q)`` in the unconstrained joint-space equations of motion
```math
M(q) \\dot{v} + c(q, v, w_\\text{ext}) = \\tau
```

This method implements the composite rigid body algorithm.
"""

"""
$(SIGNATURES)

$mass_matrix_doc

$noalloc_doc

The `out` argument must be an ``n_v \\times n_v`` lower triangular
`Symmetric` matrix, where ``n_v`` is the dimension of the `Mechanism`'s joint
velocity vector ``v``.
"""
function mass_matrix!{X, M, C}(out::Symmetric{C, Matrix{C}}, state::MechanismState{X, M, C})
    @boundscheck size(out, 1) == num_velocities(state) || error("mass matrix has wrong size")
    @boundscheck out.uplo == 'L' || error("expected a lower triangular symmetric matrix type as the mass matrix")
    fill!(out.data, zero(C))
    mechanism = state.mechanism

    for jointi in tree_joints(mechanism)
        # Hii
        irange = velocity_range(state, jointi)
        if length(irange) > 0
            bodyi = successor(jointi, mechanism)
            Si = motion_subspace_in_world(state, jointi)
            Ii = crb_inertia(state, bodyi)
            F = Ii * Si
            istart = first(irange)
            force_space_matrix_transpose_mul_jacobian!(out.data, istart, istart, F, Si, 1)

            # Hji, Hij
            body = predecessor(jointi, mechanism)
            while (!isroot(body, mechanism))
                jointj = joint_to_parent(body, mechanism)
                jrange = velocity_range(state, jointj)
                if length(jrange) > 0
                    Sj = motion_subspace_in_world(state, jointj)
                    jstart = first(jrange)
                    force_space_matrix_transpose_mul_jacobian!(out.data, istart, jstart, F, Sj, 1)
                end
                body = predecessor(jointj, mechanism)
            end
        end
    end
end

"""
$(SIGNATURES)

$mass_matrix_doc
"""
function mass_matrix{X, M, C}(state::MechanismState{X, M, C})
    nv = num_velocities(state)
    ret = Symmetric(Matrix{C}(nv, nv), :L)
    mass_matrix!(ret, state)
    ret
end

const momentum_matrix_doc = """Compute the momentum matrix ``A(q)`` of the
`Mechanism` in the given state.

The momentum matrix maps the `Mechanism`'s joint velocity vector ``v`` to
its total momentum.
"""

"""
$(SIGNATURES)

$momentum_matrix_doc

$noalloc_doc

The `out` argument must be a mutable `MomentumMatrix` with as many columns as
the dimension of the `Mechanism`'s joint velocity vector ``v``.
"""
function momentum_matrix!(out::MomentumMatrix, state::MechanismState)
    @boundscheck num_velocities(state) == num_cols(out) || error("size mismatch")
    pos = 1
    mechanism = state.mechanism
    for joint in tree_joints(mechanism)
        body = successor(joint, mechanism)
        part = crb_inertia(state, body) * motion_subspace_in_world(state, joint)
        @framecheck out.frame part.frame # TODO: transform part to out.frame instead?
        n = 3 * num_cols(part)
        copy!(out.angular, pos, part.angular, 1, n)
        copy!(out.linear, pos, part.linear, 1, n)
        pos += n
    end
end

"""
$(SIGNATURES)

$momentum_matrix_doc
"""
function momentum_matrix(state::MechanismState)
    ncols = num_velocities(state)
    T = cache_eltype(state)
    ret = MomentumMatrix(root_frame(state.mechanism), Matrix{T}(3, ncols), Matrix{T}(3, ncols))
    momentum_matrix!(ret, state)
    ret
end

function bias_accelerations!{T, X, M}(out::AbstractVector{SpatialAcceleration{T}}, state::MechanismState{X, M})
    mechanism = state.mechanism
    gravitybias = convert(SpatialAcceleration{T}, -gravitational_spatial_acceleration(mechanism))
    for joint in tree_joints(mechanism)
        body = successor(joint, mechanism)
        out[vertex_index(body)] = gravitybias + bias_acceleration(state, body)
    end
    nothing
end

function spatial_accelerations!{T, X, M}(out::AbstractVector{SpatialAcceleration{T}},
        state::MechanismState{X, M}, v̇::StridedVector)
    mechanism = state.mechanism

    # TODO: consider merging back into one loop
    # unbiased joint accelerations + gravity
    root = root_body(mechanism)
    out[vertex_index(root)] = convert(SpatialAcceleration{T}, -gravitational_spatial_acceleration(mechanism))
    for joint in tree_joints(mechanism)
        body = successor(joint, mechanism)
        S = motion_subspace_in_world(state, joint)
        v̇joint = UnsafeVectorView(v̇, velocity_range(state, joint))
        jointaccel = SpatialAcceleration(S, v̇joint)
        out[vertex_index(body)] = out[vertex_index(predecessor(joint, mechanism))] + jointaccel
    end

    # add bias acceleration - gravity
    for joint in tree_joints(mechanism)
        body = successor(joint, mechanism)
        out[vertex_index(body)] += bias_acceleration(state, body)
    end
    nothing
end

function newton_euler!{T, X, M, W}(
        out::AbstractVector{Wrench{T}}, state::MechanismState{X, M},
        accelerations::AbstractVector{SpatialAcceleration{T}},
        externalwrenches::AbstractVector{Wrench{W}})
    mechanism = state.mechanism
    for joint in tree_joints(mechanism)
        body = successor(joint, mechanism)
        wrench = newton_euler(state, body, accelerations[vertex_index(body)])
        out[vertex_index(body)] = wrench - externalwrenches[vertex_index(body)]
    end
end


function joint_wrenches_and_torques!{T, X, M}(
        torquesout::StridedVector{T},
        net_wrenches_in_joint_wrenches_out::AbstractVector{Wrench{T}},
        state::MechanismState{X, M})
    # Note: pass in net wrenches as wrenches argument. wrenches argument is modified to be joint wrenches
    @boundscheck length(torquesout) == num_velocities(state) || error("length of torque vector is wrong")
    mechanism = state.mechanism
    joints = tree_joints(mechanism)
    for i = length(joints) : -1 : 1
        joint = joints[i]
        body = successor(joint, mechanism)
        parentbody = predecessor(joint, mechanism)
        jointwrench = net_wrenches_in_joint_wrenches_out[vertex_index(body)]
        if !isroot(parentbody, mechanism)
            # TODO: consider also doing this for the root:
            net_wrenches_in_joint_wrenches_out[vertex_index(parentbody)] += jointwrench # action = -reaction
        end
        jointwrench = transform(jointwrench, inv(transform_to_root(state, body))) # TODO: stay in world frame?
        @inbounds τjoint = UnsafeVectorView(torquesout, velocity_range(state, joint))
        joint_torque!(joint, τjoint, configuration(state, joint), jointwrench)
    end
end

"""
$(SIGNATURES)

Compute the 'dynamics bias term', i.e. the term
```math
c(q, v, w_\\text{ext})
```
in the unconstrained joint-space equations of motion
```math
M(q) \\dot{v} + c(q, v, w_\\text{ext}) = \\tau
```

$noalloc_doc
"""
function dynamics_bias!{T, X, M, W}(
        torques::AbstractVector{T},
        biasaccelerations::AbstractVector{SpatialAcceleration{T}},
        wrenches::AbstractVector{Wrench{T}},
        state::MechanismState{X, M},
        externalwrenches::AbstractVector{Wrench{W}} = ConstVector(zero(Wrench{T}, root_frame(state.mechanism)), num_bodies(state.mechanism)))
    bias_accelerations!(biasaccelerations, state)
    newton_euler!(wrenches, state, biasaccelerations, externalwrenches)
    joint_wrenches_and_torques!(torques, wrenches, state)
end

const inverse_dynamics_doc = """Do inverse dynamics, i.e. compute ``\\tau``
in the unconstrained joint-space equations of motion
```math
M(q) \\dot{v} + c(q, v, w_\\text{ext}) = \\tau
```
given joint configuration vector ``q``, joint velocity vector ``v``,
joint acceleration vector ``\\dot{v}`` and (optionally) external
wrenches ``w_\\text{ext}``.

This method implements the recursive Newton-Euler algorithm.

Currently doesn't support `Mechanism`s with cycles.
"""

"""
$(SIGNATURES)

$inverse_dynamics_doc

$noalloc_doc
"""
function inverse_dynamics!{T, X, M, V, W}(
        torquesout::AbstractVector{T},
        jointwrenchesout::AbstractVector{Wrench{T}},
        accelerations::AbstractVector{SpatialAcceleration{T}},
        state::MechanismState{X, M},
        v̇::AbstractVector{V},
        externalwrenches::AbstractVector{Wrench{W}} = ConstVector(zero(Wrench{T}, root_frame(state.mechanism)), num_bodies(state.mechanism)))
    length(tree_joints(state.mechanism)) == length(joints(state.mechanism)) || error("This method can currently only handle tree Mechanisms.")
    spatial_accelerations!(accelerations, state, v̇)
    newton_euler!(jointwrenchesout, state, accelerations, externalwrenches)
    joint_wrenches_and_torques!(torquesout, jointwrenchesout, state)
end

"""
$(SIGNATURES)

$inverse_dynamics_doc
"""
function inverse_dynamics{X, M, V, W}(
        state::MechanismState{X, M},
        v̇::AbstractVector{V},
        externalwrenches::AbstractVector{Wrench{W}} = ConstVector(zero(Wrench{X}, root_frame(state.mechanism)), num_bodies(state.mechanism)))
    T = promote_type(X, M, V, W)
    torques = Vector{T}(num_velocities(state))
    nb = num_bodies(state.mechanism)
    jointwrenches = Vector{Wrench{T}}(nb)
    accelerations = Vector{SpatialAcceleration{T}}(nb)
    inverse_dynamics!(torques, jointwrenches, accelerations, state, v̇, externalwrenches)
    torques
end

function constraint_jacobian_and_bias!(state::MechanismState, constraintjacobian::AbstractMatrix, constraintbias::AbstractVector)
    # TODO: allocations
    structure = state.constraint_jacobian_structure
    mechanism = state.mechanism
    rows = rowvals(structure)
    signs = nonzeros(structure)
    nontreejoints = non_tree_joints(state)
    rowstart = 1
    constraintjacobian[:] = 0
    for i = 1 : length(nontreejoints)
        nontreejoint = nontreejoints[i]
        nextrowstart = rowstart + num_constraints(nontreejoint)
        range = rowstart : nextrowstart - 1

        # Constraint wrench subspace.
        jointtransform = relative_transform(state, frame_after(nontreejoint), frame_before(nontreejoint)) # TODO: expensive
        T = constraint_wrench_subspace(nontreejoint, jointtransform)
        T = transform(T, transform_to_root(state, T.frame)) # TODO: expensive

        # Jacobian rows.
        for j in nzrange(structure, i)
            treejoint = tree_joints(mechanism)[rows[j]]
            J = motion_subspace_in_world(state, treejoint)
            sign = signs[j]
            colstart = first(velocity_range(state, treejoint))
            force_space_matrix_transpose_mul_jacobian!(constraintjacobian, rowstart, colstart, T, J, sign)
        end

        # Constraint bias.
        has_fixed_subspaces(nontreejoint) || error("Only joints with fixed motion subspace (Ṡ = 0) supported at this point.")
        kjoint = UnsafeVectorView(constraintbias, range)
        pred = predecessor(nontreejoint, mechanism)
        succ = successor(nontreejoint, mechanism)
        biasaccel = cross(twist_wrt_world(state, succ), twist_wrt_world(state, pred)) + (bias_acceleration(state, succ) + -bias_acceleration(state, pred)) # 8.47 in Featherstone
        At_mul_B!(kjoint, T, biasaccel)
        rowstart = nextrowstart
    end
end

function dynamics_solve!(result::DynamicsResult, τ::AbstractVector)
    # version for general scalar types
    # TODO: make more efficient
    M = result.massmatrix
    c = result.dynamicsbias
    v̇ = result.v̇

    K = result.constraintjacobian
    k = result.constraintbias
    λ = result.λ

    nv = size(M, 1)
    nl = size(K, 1)
    G = [M K';
         K zeros(nl, nl)]
    r = [τ - c; -k]
    v̇λ = G \ r
    v̇[:] = view(v̇λ, 1 : nv)
    λ[:] = view(v̇λ, nv + 1 : nv + nl)
    nothing
end

function dynamics_solve!{S<:Number, T<:LinAlg.BlasReal}(result::DynamicsResult{S, T}, τ::AbstractVector{T})
    # optimized version for BLAS floats
    M = result.massmatrix
    c = result.dynamicsbias
    v̇ = result.v̇

    K = result.constraintjacobian
    k = result.constraintbias
    λ = result.λ

    L = result.L
    A = result.A
    z = result.z
    Y = result.Y

    L[:] = M.data
    uplo = M.uplo
    LinAlg.LAPACK.potrf!(uplo, L) # L <- Cholesky decomposition of M; M == L Lᵀ (note: Featherstone, page 151 uses M == Lᵀ L instead)
    τBiased = v̇
    @simd for i in eachindex(τ)
        @inbounds τBiased[i] = τ[i] - c[i]
    end
    # TODO: use τBiased .= τ .- c in 0.6

    if size(K, 1) > 0
        # Loops.
        # Basic math:
        # DAE:
        #   [M Kᵀ] [v̇] = [τ - c]
        #   [K 0 ] [λ]  = [-k]
        # Solve for v̇:
        #   v̇ = M⁻¹ (τ - c - Kᵀ λ)
        # Plug into λ equation:
        #   K M⁻¹ (τ - c - Kᵀ λ) = -k
        #   K M⁻¹ Kᵀ λ = K M⁻¹(τ - c) + k
        # which can be solved for λ, after which λ can be used to solve for v̇.

        # Implementation loosely follows page 151 of Featherstone, Rigid Body Dynamics Algorithms.
        # It is slightly different because Featherstone uses M = Lᵀ L, whereas BLAS uses M = L Lᵀ.
        # M = L Lᵀ (Cholesky decomposition)
        # Y = K L⁻ᵀ, so that Y Yᵀ = K L⁻ᵀ L⁻¹ Kᵀ = K M⁻¹ Kᵀ = A
        # z = L⁻¹ (τ - c)
        # b = Y z + k
        # solve A λ = b for λ
        # solve M v̇  = τ - c for v̇

        # Compute Y = K L⁻ᵀ
        Y[:] = K
        LinAlg.BLAS.trsm!('R', uplo, 'T', 'N', one(T), L, Y)

        # Compute z = L⁻¹ (τ - c)
        z[:] = τBiased
        LinAlg.BLAS.trsv!(uplo, 'N', 'N', L, z) # z <- L⁻¹ (τ - c)

        # Compute A = Y Yᵀ == K * M⁻¹ * Kᵀ
        LinAlg.BLAS.gemm!('N', 'T', one(T), Y, Y, zero(T), A) # A <- K * M⁻¹ * Kᵀ

        # Compute b = Y z + k
        b = λ
        b[:] = k
        LinAlg.BLAS.gemv!('N', one(T), Y, z, one(T), b) # b <- Y z + k

        # Compute λ = A⁻¹ b == (K * M⁻¹ * Kᵀ)⁻¹ * (K * M⁻¹ * (τ - c) + k)
        LinAlg.LAPACK.posv!(uplo, A, b) # b == λ <- (K * M⁻¹ * Kᵀ)⁻¹ * (K * M⁻¹ * (τ - c) + k)

        # Update τBiased: subtract Kᵀ * λ
        LinAlg.BLAS.gemv!('T', -one(T), K, λ, one(T), τBiased) # τBiased <- τ - c - Kᵀ * λ

        # Solve for v̇ = M⁻¹ * (τ - c - Kᵀ * λ)
        LinAlg.LAPACK.potrs!(uplo, L, τBiased) # τBiased ==v̇ <- M⁻¹ * (τ - c - Kᵀ * λ)
    else
        # No loops.
        LinAlg.LAPACK.potrs!(uplo, L, τBiased) # τBiased == v̇ <- M⁻¹ * (τ - c)
    end
    nothing
end

"""
$(SIGNATURES)

Compute the joint acceleration vector ``\\dot{v}`` and Lagrange multipliers
``\\lambda`` that satisfy the joint-space equations of motion
```math
M(q) \\dot{v} + c(q, v, w_\\text{ext}) = \\tau - K(q)^{T} \\lambda
```
and the constraint equations
```math
K(q) \\dot{v} = -k
```
given joint configuration vector ``q``, joint velocity vector ``v``, and
(optionally) joint torques ``\\tau`` and external wrenches ``w_\\text{ext}``.
"""
function dynamics!{T, X, M, Tau, W}(result::DynamicsResult{T}, state::MechanismState{X, M},
        torques::AbstractVector{Tau} = ConstVector(zero(T), num_velocities(state)),
        externalwrenches::AbstractVector{Wrench{W}} = ConstVector(zero(Wrench{T}, root_frame(state.mechanism)), num_bodies(state.mechanism)))
    fill!(result.ṡ, zero(T)) # TODO
    dynamics_bias!(result.dynamicsbias, result.accelerations, result.jointwrenches, state, externalwrenches)
    mass_matrix!(result.massmatrix, state)
    constraint_jacobian_and_bias!(state, result.constraintjacobian, result.constraintbias)
    dynamics_solve!(result, torques)
    nothing
end


"""
$(SIGNATURES)

Convenience function for use with standard ODE integrators that takes a
`Vector` argument
```math
x = \\left(\\begin{array}{c}
q\\\\
v
\\end{array}\\right)
```
and returns a `Vector` ``\\dot{x}``.
"""
function dynamics!{T, X, M, Tau, W}(ẋ::StridedVector{X},
        result::DynamicsResult{T}, state::MechanismState{X, M}, stateVec::AbstractVector{X},
        torques::AbstractVector{Tau} = ConstVector(zero(T), num_velocities(state)),
        externalwrenches::AbstractVector{Wrench{W}} = ConstVector(zero(Wrench{T}, root_frame(state.mechanism)), num_bodies(state.mechanism)))
    set!(state, stateVec)
    nq = num_positions(state)
    nv = num_velocities(state)
    q̇ = view(ẋ, 1 : nq) # allocates
    v̇ = view(ẋ, nq + 1 : nq + nv) # allocates
    configuration_derivative!(q̇, state)
    dynamics!(result, state, torques, externalwrenches)
    copy!(v̇, result.v̇)
    ẋ
end
