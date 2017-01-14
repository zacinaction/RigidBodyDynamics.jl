module Segfault
using StaticArrays

typealias ContiguousSMatrixColumnView{S1, S2, T, L} SubArray{T,2,SMatrix{S1, S2, T, L},Tuple{Colon,UnitRange{Int64}},true}


function smatrix3x6view{T}(mat::SMatrix{3, 0, T, 0})
    data = fill(NaN, SMatrix{3, 6, T, 18})
    ret = view(data, :, 1 : 0)
end

immutable WrenchMatrix{A<:AbstractMatrix}
    angular::A
    linear::A

    function WrenchMatrix(angular::A, linear::A)
        new(angular, linear)
    end
end

function WrenchMatrix{A<:AbstractMatrix}(angular::A, linear::A)
    WrenchMatrix{A}(angular, linear)
end

typealias WrenchSubspace{T} WrenchMatrix{ContiguousSMatrixColumnView{3, 6, T, 18}}

abstract JointType{T<:Real}
immutable QuaternionFloating{T} <: JointType{T}
end
function _constraint_wrench_subspace{T}(::QuaternionFloating{T})
    angular = zeros(SMatrix{3, 0, T})
    linear = zeros(SMatrix{3, 0, T})
    WrenchMatrix(smatrix3x6view(angular), smatrix3x6view(linear))
end

type Joint{T<:Real}
    jointType::JointType{T}
end

function constraint_wrench_subspace{M}(joint::Joint{M})::WrenchSubspace{M}
    _constraint_wrench_subspace(joint.jointType::QuaternionFloating{M})
end

end # module

js = (Segfault.Joint(Segfault.QuaternionFloating{Float64}()),)
for joint in js
    T = Segfault.constraint_wrench_subspace(joint)
    println(length(T.angular))
end
