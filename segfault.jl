module Segfault
using StaticArrays

# typealias for typeof(view(rand(SMatrix{S1, S2, T, L}), :, 1:0)):
typealias ViewType{S1, S2, T, L} SubArray{T,2,SMatrix{S1, S2, T, L},Tuple{Colon,UnitRange{Int64}},true} # for 0.5
# typealias ViewType{S1, S2, T, L} SubArray{T,2,SMatrix{S1,S2,T,L},Tuple{Base.Slice{Base.OneTo{Int64}},UnitRange{Int64}},true} # for latest master

immutable MatrixHolder{A<:AbstractMatrix}
    a::A
end

typealias SpecificMatrixHolder{T} MatrixHolder{ViewType{3, 6, T, 18}}

type Foo{T} end

function bar{T}(fooHolder::Foo{T})::SpecificMatrixHolder{T} # return type annotation required
    data = fill(NaN, SMatrix{3, 6, T})
    MatrixHolder(view(data, :, 1 : 0))
end

end # module

fooHolders = (Segfault.Foo{Float64}(),)
for fooHolder in fooHolders
    arrayHolder = Segfault.bar(fooHolder) # call to bar must occur in for loop
    println(length(arrayHolder.a))
end
