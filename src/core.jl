# These functions are taken from algorithms.jl in
# https://github.com/JuliaImages/Images.jl/. We don't want to depend on the
# entire Images package just for these functions hence I have just copied and
# pasted them here for now. However, as we are refactoring the Images.jl package
# we might want to move these (and some other functions) into a separate
# package. Once that is completed we can import that package and remove this
# file.

"""
`m = minfinite(A)` calculates the minimum value in `A`, ignoring any values that are not finite (Inf or NaN).
"""
function minfinite(A::AbstractArray{T}) where T
    ret = sentinel_min(T)
    for a in A
        ret = minfinite_scalar(a, ret)
    end
    ret
end
function minfinite(f, A::AbstractArray)
    ret = sentinel_min(typeof(f(first(A))))
    for a in A
        ret = minfinite_scalar(f(a), ret)
    end
    ret
end

"""
`m = maxfinite(A)` calculates the maximum value in `A`, ignoring any values that are not finite (Inf or NaN).
"""
function maxfinite(A::AbstractArray{T}) where T
    ret = sentinel_max(T)
    for a in A
        ret = maxfinite_scalar(a, ret)
    end
    ret
end
function maxfinite(f, A::AbstractArray)
    ret = sentinel_max(typeof(f(first(A))))
    for a in A
        ret = maxfinite_scalar(f(a), ret)
    end
    ret
end

minfinite_scalar(a::T, b::T) where {T} = isfinite(a) ? (b < a ? b : a) : b
maxfinite_scalar(a::T, b::T) where {T} = isfinite(a) ? (b > a ? b : a) : b
minfinite_scalar(a::T, b::T) where {T<:Union{Integer,FixedPoint}} = b < a ? b : a
maxfinite_scalar(a::T, b::T) where {T<:Union{Integer,FixedPoint}} = b > a ? b : a
minfinite_scalar(a, b) = minfinite_scalar(promote(a, b)...)
maxfinite_scalar(a, b) = maxfinite_scalar(promote(a, b)...)

function minfinite_scalar(c1::C, c2::C) where C<:AbstractRGB
    C(minfinite_scalar(c1.r, c2.r),
      minfinite_scalar(c1.g, c2.g),
      minfinite_scalar(c1.b, c2.b))
end
function maxfinite_scalar(c1::C, c2::C) where C<:AbstractRGB
    C(maxfinite_scalar(c1.r, c2.r),
      maxfinite_scalar(c1.g, c2.g),
      maxfinite_scalar(c1.b, c2.b))
end

sentinel_min(::Type{T}) where {T<:Union{Integer,FixedPoint}} = typemax(T)
sentinel_max(::Type{T}) where {T<:Union{Integer,FixedPoint}} = typemin(T)
sentinel_min(::Type{T}) where {T<:AbstractFloat} = convert(T, NaN)
sentinel_max(::Type{T}) where {T<:AbstractFloat} = convert(T, NaN)
sentinel_min(::Type{C}) where {C<:AbstractRGB} = _sentinel_min(C, eltype(C))
_sentinel_min(::Type{C},::Type{T}) where {C<:AbstractRGB,T} = (s = sentinel_min(T); C(s,s,s))
sentinel_max(::Type{C}) where {C<:AbstractRGB} = _sentinel_max(C, eltype(C))
_sentinel_max(::Type{C},::Type{T}) where {C<:AbstractRGB,T} = (s = sentinel_max(T); C(s,s,s))
sentinel_min(::Type{C}) where {C<:AbstractGray} = _sentinel_min(C, eltype(C))
_sentinel_min(::Type{C},::Type{T}) where {C<:AbstractGray,T} = C(sentinel_min(T))
sentinel_max(::Type{C}) where {C<:AbstractGray} = _sentinel_max(C, eltype(C))
_sentinel_max(::Type{C},::Type{T}) where {C<:AbstractGray,T} = C(sentinel_max(T))
