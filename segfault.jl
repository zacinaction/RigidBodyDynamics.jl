module Segfault
using StaticArrays

# typealias for typeof(view(rand(SMatrix{S1, S2, T, L}), :, 1:0)):
typealias ViewType{S1, S2, T, L} SubArray{T,2,SMatrix{S1, S2, T, L},Tuple{Colon,UnitRange{Int64}},true}

immutable ArrayHolder{A<:AbstractMatrix}
    a::A
end

typealias SpecificArrayHolder{T} ArrayHolder{ViewType{3, 6, T, 18}}

immutable Foo{T} end

type FooHolder{T<:Real}
    foo::Foo{T}
end

function smatrix3x6view(mat::SMatrix{3, 1}) # mat's type must not be concrete
    T = eltype(mat)
    data = fill(NaN, SMatrix{3, 6, T})
    ret = view(data, :, 1 : 0)
end

function _bar{T}(::Foo{T})
    a = zeros(SMatrix{3, 1, T})
    ArrayHolder(smatrix3x6view(a))
end

function bar{M}(fooHolder::FooHolder{M})::SpecificArrayHolder{M} # return type annotation required
    _bar(fooHolder.foo)
end

end # module

fooHolders = (Segfault.FooHolder(Segfault.Foo{Float64}()),)
for fooHolder in fooHolders
    arrayHolder = Segfault.bar(fooHolder) # call to bar must occur in for loop
    println(length(arrayHolder.a))
end
