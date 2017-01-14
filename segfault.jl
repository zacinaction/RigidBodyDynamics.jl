module Segfault
using StaticArrays, Rotations

import Base: rand, +

typealias ContiguousSMatrixColumnView{S1, S2, T, L} SubArray{T,2,SMatrix{S1, S2, T, L},Tuple{Colon,UnitRange{Int64}},true}


@generated function smatrix3x6view{N, T}(mat::SMatrix{3, N, T})
    colrange = 1 : N
    if N == 0
        data = fill(NaN, SMatrix{3, 6, T, 18})
        ret = view(data, :, 1 : 0)
        return :($ret)
    elseif N == 6
        return :(ContiguousSMatrixColumnView{3, 6, $T, 18}(mat, (:, $colrange), 0, 1)) # faster way to create view)
    else
        fillerSize = 6 - N
        filler = fill(NaN, SMatrix{3, fillerSize, T, 3 * fillerSize})
        return quote
            data = hcat(mat, $filler)::SMatrix{3, 6, $T, 18}
            ContiguousSMatrixColumnView{3, 6, $T, 18}(data, (:, $colrange), 0, 1)
        end
    end
end

macro rtti_dispatch(typeTuple, signature)
    @assert signature.head == :call
    @assert length(signature.args) > 1
    @assert typeTuple.head == :tuple

    f = signature.args[1]
    args = signature.args[2 : end]
    dispatchArg = args[1]
    otherArgs = args[2 : end]
    types = typeTuple.args

    ret = :(error("type not recognized"))
    for T in reverse(types)
        ret = Expr(:if, :(isa($dispatchArg, $T)), :(return $(f)($(dispatchArg)::$T, $(otherArgs...))), ret)
    end
    :($(esc(ret)))
end

const next_frame_id = Ref(0)
const frame_names = Dict{Int64, String}()

immutable CartesianFrame3D
    id::Int64

    function CartesianFrame3D(name::String)
        ret = new(next_frame_id.x)
        next_frame_id.x += 1
        frame_names[ret.id] = name
        ret
    end

    function CartesianFrame3D()
        ret = new(next_frame_id.x)
        next_frame_id.x += 1
        ret
    end
end

immutable WrenchMatrix{A<:AbstractMatrix}
    frame::CartesianFrame3D
    angular::A
    linear::A

    function WrenchMatrix(frame::CartesianFrame3D, angular::A, linear::A)
        @assert size(angular, 1) == 3
        @assert size(linear, 1) == 3
        @assert size(angular, 2) == size(linear, 2)
        new(frame, angular, linear)
    end
end

function WrenchMatrix{A<:AbstractMatrix}(frame::CartesianFrame3D, angular::A, linear::A)
    WrenchMatrix{A}(frame, angular, linear)
end

typealias WrenchSubspace{T} WrenchMatrix{ContiguousSMatrixColumnView{3, 6, T, 18}}
function WrenchSubspace(frame::CartesianFrame3D, angular, linear)
    WrenchMatrix(frame, smatrix3x6view(angular), smatrix3x6view(linear))
end

abstract JointType{T<:Real}
immutable QuaternionFloating{T} <: JointType{T}
end
function _constraint_wrench_subspace{T<:Real, X<:Real}(::QuaternionFloating{T}, frameAfter::CartesianFrame3D, q::AbstractVector{X})
    S = promote_type(T, X)
    WrenchSubspace(frameAfter, zeros(SMatrix{3, 0, S}), zeros(SMatrix{3, 0, S}))
end

type Joint{T<:Real}
    frameAfter::CartesianFrame3D
    jointType::JointType{T}

    Joint(name::String, jointType::JointType{T}) = new(CartesianFrame3D(string("after_", name)), jointType)
end

Joint{T<:Real}(name::String, jointType::JointType{T}) = Joint{T}(name, jointType)


function constraint_wrench_subspace{M, X}(joint::Joint{M}, q::AbstractVector{X})::WrenchSubspace{promote_type(M, X)}
    @rtti_dispatch (QuaternionFloating{M},) _constraint_wrench_subspace(joint.jointType, joint.frameAfter, q)
end

end # module

js = (Segfault.Joint("bla", Segfault.QuaternionFloating{Float64}()),)
for joint in js
    T = Segfault.constraint_wrench_subspace(joint, rand(7))
    println(length(T.angular))
end
