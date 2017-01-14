module Segfault
using StaticArrays

typealias ContiguousSMatrixColumnView{S1, S2, T, L} SubArray{T,2,SMatrix{S1, S2, T, L},Tuple{Colon,UnitRange{Int64}},true}

immutable WrenchMatrix{A<:AbstractMatrix}
    angular::A
end

typealias WrenchSubspace{T} WrenchMatrix{ContiguousSMatrixColumnView{3, 6, T, 18}}

abstract JointType{T<:Real}

immutable QuaternionFloating{T} <: JointType{T} end

type Joint{T<:Real}
    jointType::JointType{T}
end

function smatrix3x6view(mat::SMatrix{3, 0}) # mat must not have concrete type
    T = eltype(mat)
    data = fill(NaN, SMatrix{3, 6, T})
    ret = view(data, :, 1 : 0)
end

function _constraint_wrench_subspace{T}(::QuaternionFloating{T})
    angular = zeros(SMatrix{3, 0, T})
    WrenchMatrix(smatrix3x6view(angular))
end

function constraint_wrench_subspace{M}(joint::Joint{M})::WrenchSubspace{M}
    _constraint_wrench_subspace(joint.jointType::QuaternionFloating{M})
end

end # module

joints = (Segfault.Joint(Segfault.QuaternionFloating{Float64}()),)
for joint in joints
    T = Segfault.constraint_wrench_subspace(joint)
    println(length(T.angular))
end
