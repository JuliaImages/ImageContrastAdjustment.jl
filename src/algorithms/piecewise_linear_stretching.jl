"""
```
    PiecewiseLinearStretching <: AbstractHistogramAdjustmentAlgorithm
    PiecewiseLinearStretching((0.0, 0.2, 0.8, 1.0), (0.0, 0.1, 0.7, 1.0) ; saturate = true)
    PiecewiseLinearStretching([0.0, 0.2, 0.8, 1.0], [0.0, 0.1, 0.7, 1.0] ; saturate = true)
    PiecewiseLinearStretching(; src_knots = (0.0, 0.2, 0.8, 1.0), 
                                dst_knots = (0.0, 0.1, 0.7, 1.0), saturate = true)
    PiecewiseLinearStretching(ClosedInterval(0.0, 0.2) => ClosedInterval(0.0, 0.1),
                              Interval{:open, :closed}(0.2, 0.8) => Interval{:open, :closed}(0.1, 0.7),
                              Interval{:open, :closed}(0.8, 1,0) => Interval{:open, :closed}(0.7, 1.0))
    PiecewiseLinearStretching(Percentiles(1,99), (0,1), img::GenericGrayImage)
    PiecewiseLinearStretching(Percentiles(1,50,99), (0, 0.5, 1), img::GenericGrayImage)
    PiecewiseLinearStretching(MinMax(), (0,1), img::GenericGrayImage)                          
    adjust_histogram([T,] img, f::PiecewiseLinearStretching)
    adjust_histogram!([out,] img, f::PiecewiseLinearStretching)
```

When used with `adjust_histogram`, returns an image for which the intensity range was
partitioned into subintervals  and mapped linearly to a new set of subintervals. 

# Details

Piecewise linear stretching is a contrast enhancing transformation that is used to modify
the dynamic range of an image. The concept involves partitioning the intensity
range of an image, ``I = [A,B]`` into a sequence of ``m`` disjoint intervals

```math
I_1 = [A_1,A_2], I_2 = (A_2, A_3], \\ldots, I_{m} = (A_m,A_{m+1}]
```
which encompass the entire range of intensity values such that ``\\sum_j |I_j | = | I |``. 
One subsequently constructs a concomitant new set of disjoint intervals
```math
I_{1}^{′} = [a_1,a_2], I_{2}^{′} = (a_2, a_3], \\ldots, I_{m}^{′} = (a_m,a_{m+1}]
```
such that ``\\sum_j |I_{j}^{′} | = | I^{′} |``, where ``I^{′} = [a,b]``. Finally, one introduces a
sequence of mappings ``f_j : I_j \\to I_{j}^{′}`` ``(j = 1, \\ldots, m)`` given by
```math
f_j(x) = (x-A_{j}) \\frac{a_{j+1}-a_j}{A_{j+1}-A_{j}} + a_j.
```
which linearly map intensities from the interval ``I_j`` to the interval ``I_{j}^{′}`` .

# Options

Various options for the parameters of the `adjust_histogram` and
`PiecewiseLinearStretching` type are described in more detail below.

## Choices for `img`

The function can handle a variety of input types. The returned image depends on the input
type.

For colored images, the input is converted to the [YIQ](https://en.wikipedia.org/wiki/YIQ)
type and the intensities of the Y channel are stretched to the specified range. The modified
Y channel is then combined with the I and Q channels and the resulting image converted to
the same type as the input.

## Specifying intervals explicitly
One can specify the mappings between source and destination intervals using
the interval types provided by [IntervalSets.jl](https://github.com/JuliaMath/IntervalSets.jl).
For instance, the mappings [0, 0.2] => [0, 0.1], (0.2, 0.8] => (0.1, 0.7] and (0.8, 1.0] =>
(0.7, 1.0] can be specified as follows:
```julia
PiecewiseLinearStretching(ClosedInterval(0.0, 0.2) => ClosedInterval(0.0, 0.1),
                          Interval{:open, :closed}(0.2, 0.8) => Interval{:open, :closed}(0.1, 0.7),
                          Interval{:open, :closed}(0.8, 1,0) => Interval{:open, :closed}(0.7, 1.0))
```
While this manner of specifying intervals permits the most flexibility, it is
cumbersome for many types of mappings one would like to specify in practice. 
Hence, one can instead specify the mappings by designating a sequence of knot pairs. 

## Specifying intervals implicitly as knots
As an alternative to specifying the endpoints of each interval, one can list the knots of
the intervals. For instance, the mappings [0.0, 0.2] => [0, 0.1], (0.2, 0.8] => (0.1, 0.7]
and (0.8, 1.0] => (0.7, 1.0] can be listed as follows:

```julia
# Without the use of keyword arguments.
PiecewiseLinearStretching((0.0, 0.2, 0.8, 1.0), (0.0, 0.1, 0.7, 1.0))
# With the use of keyword arguments.
PiecewiseLinearStretching(src_knots = (0.0, 0.2, 0.8, 1.0), dst_knots = (0.0, 0.1, 0.7, 1.0))
```

## Specifying intervals implicitly as [`Percentiles`](@ref).
Instead of directly listing knots, one can set them implicitly as percentiles of the 
image pixels. For example, to may intensities between the 1st and 99th percentile to the
unit interval one could use the following code:
```julia
img = testimage("mandril_gray")
# Note that you have to pass in the img as an argument to obtain the percentiles. 
f = PiecewiseLinearStretching(Percentile((1,99)), (0,1), img) 
img_adjusted = adjust_histogram(img, f)
```

## Specifing intervals implicitly as [`MinMax`](@ref).
One can also specify and interval implicitly as the minimum and maximum intensity in an
image. For example, one can map intensities within the minimum/maximum intensity in an image
to the unit interval:
```julia
img = testimage("mandril_gray")
# Note that you have to pass in the img as an argument to obtain the minimum and maximum intensities. 
f = PiecewiseLinearStretching(MinMax(), (0,1), img) 
img_adjusted = adjust_histogram(img, f)
```

## Constraints on intervals and knots. 
To specify a single interval you must include a trailing comma, e.g.
```julia
# Using the ellipses notation provided by IntervalSets.jl
f = PiecewiseLinearStretching((0.0..1.0 => 0.2..1.0,))
# Using the types provided by IntervalSets.jl
f = PiecewiseLinearStretching(ClosedInterval(0.0, 1.0) => ClosedInterval(0.2, 1.0),)
```

The endpoints of any source interval cannot be equal. For example, `ClosedInterval(0.5, 0.5)
=> ClosedInterval(0.0, 1.0)` is not permitted because the value `0.5` cannot be
simultaneously mapped to `0.0` and `1.0`. There is no such constraint for the destination
interval. For instance, `ClosedInterval(0.0, 1.0) => ClosedInterval(0.5, 0.5)` is permitted
and corresponds to mapping all intensities in the unit interval to the value `0.5`. 

Care must be taken to ensure that the destination subintervals  are representable within the
type of image on which the piecewise linear stretching will be applied. For example, if you
specify a negative output interval like `(0.0..1.0 => -1.0..0,),` and try to apply the
transformation on an image with type `Gray{N0f8}`, you will receive a DomainError. 

When listing knots, the length of `src_knots` must match the length of `dst_knots`.


## Saturating pixels

When listing knots, one has the option of saturating pixels that fall outside the endpoints
of the source knots to the endpoints of the destination knots (this option is enabled by default). For example,
```julia
f = PiecewiseLinearStretching(src_knots = (0.2, 0.8), dst_knots = (0.0, 1.0); saturate = true)
```
will ensure that any image intensities below `0.2` are set to zero, and any intensities above
`0.8` are set to one, whilst also linearly stretching the interval `[0.2, 0.8]` to `[0.0. 1.0]`.


# Example

```julia
using ImageContrastAdjustment, TestImages

img = testimage("mandril_gray")
f = PiecewiseLinearStretching((0.0, 0.2, 0.8, 1.0), (0.0, 0.1, 0.7, 1.0))
imgo = adjust_histogram(img, f)
```

# References
1. W. Burger and M. J. Burge. *Digital Image Processing*. Texts in Computer Science, 2016. [doi:10.1007/978-1-4471-6684-9](https://doi.org/10.1007/978-1-4471-6684-9)
"""
struct PiecewiseLinearStretching{T} <: AbstractHistogramAdjustmentAlgorithm where T <: Tuple{Vararg{Pair}}
    intervals::T 
    function PiecewiseLinearStretching(intervals::T₁) where T₁ <: Tuple{Vararg{Pair}} 
        for (src, _) in intervals
            if leftendpoint(src) == rightendpoint(src) 
                throw(ArgumentError("The start and end of a src interval cannot be equal."))
            end
        end
        new{T₁}(intervals)
    end    
end

function PiecewiseLinearStretching(; src_knots, dst_knots, saturate::Bool = true)
   return PiecewiseLinearStretching(src_knots, dst_knots; saturate = saturate)
end

function PiecewiseLinearStretching(img::GenericGrayImage; src_knots, dst_knots, saturate::Bool = true)
    return PiecewiseLinearStretching(src_knots, dst_knots, img; saturate = saturate)
end

function PiecewiseLinearStretching(src_knots::MinMax, dst_knots::T, img::GenericGrayImage; saturate::Bool = true) where T <: Union{Tuple{Vararg{Union{Real,AbstractGray}}}, AbstractVector, Percentiles}
    img_min, img_max = minfinite(img), maxfinite(img)  
    if (img_min ≈ img_max)
        throw(DomainError("The image consists of a single intensity, which gives rise to a degenerate source interval."))
    end  
    if T <: Percentiles
        return PiecewiseLinearStretching((img_min, img_max), dst_knots, img; saturate = saturate)
    else
        return PiecewiseLinearStretching((img_min, img_max), dst_knots; saturate = saturate)
    end
end

function PiecewiseLinearStretching(src_knots::T, dst_knots::MinMax, img::GenericGrayImage; saturate::Bool = true) where T <: Union{Tuple{Vararg{Union{Real,AbstractGray}}}, AbstractVector, Percentiles}
    img_min, img_max = minfinite(img), maxfinite(img)  
    if T <: Percentiles
        return PiecewiseLinearStretching(src_knots, (img_min, img_max), img ; saturate = saturate)
    else
        return PiecewiseLinearStretching(src_knots, (img_min, img_max); saturate = saturate)
    end
end

function PiecewiseLinearStretching(src_knots::AbstractVector, dst_knots::Tuple{Vararg{Union{Real,AbstractGray}}}; saturate::Bool = true)
    return PiecewiseLinearStretching(tuple(src_knots...), dst_knots; satureate = saturate)
end

function PiecewiseLinearStretching(src_knots::Tuple{Vararg{Union{Real,AbstractGray}}}, dst_knots::AbstractVector; saturate::Bool = true)
    return PiecewiseLinearStretching(src_knots, tuple(dst_knots...); saturate = saturate)
end

function PiecewiseLinearStretching(src_percentile::Percentiles, dst_percentile::Percentiles, img::GenericGrayImage; saturate::Bool = true) 
    if length(src_percentile.p) != length(dst_percentile.p)
        throw(ArgumentError("The number of src percentiles and dst percentiles must match."))
    end
    src_knots = percentile(vec(img), collect(src_percentile.p))
    dst_knots = percentile(vec(img), collect(dst_percentile.p))
    return PiecewiseLinearStretching(tuple(src_knots...), tuple(dst_knots...); saturate = saturate)
end

function PiecewiseLinearStretching(src_percentile::Percentiles, dst_knots::T, img::GenericGrayImage; saturate::Bool = true) where T <: Union{Tuple{Vararg{Union{Real,AbstractGray}}}, AbstractVector}
    if length(src_percentile.p) != length(dst_knots)
        throw(ArgumentError("The number of src percentiles and destination knots must match."))
    end
    src_knots = percentile(vec(img), collect(src_percentile.p))
    return PiecewiseLinearStretching(tuple(src_knots...), dst_knots; saturate = saturate)
end

function PiecewiseLinearStretching(src_knots::T, dst_percentile::Percentiles, img::GenericGrayImage; saturate::Bool = true) where T <: Union{Tuple{Vararg{Union{Real,AbstractGray}}}, AbstractVector}
    if length(src_knots) != length(dst_percentile.p)
        throw(ArgumentError("The number of src knots and destination percentiles must match."))
    end
    dst_knots = percentile(vec(img), collect(dst_percentile.p))
    return PiecewiseLinearStretching(src_knots, tuple(dst_knots...); saturate = saturate)
end

function PiecewiseLinearStretching(src_knots::AbstractVector, dst_knots::AbstractVector; saturate::Bool = true)
    return PiecewiseLinearStretching(tuple(src_knots...), tuple(dst_knots...); saturate = saturate)
end

function PiecewiseLinearStretching(src_knots::Tuple{Vararg{Union{Real,AbstractGray}}}, dst_knots::Tuple{Vararg{Union{Real,AbstractGray}}}; saturate::Bool = true) 
    if length(src_knots) != length(dst_knots)
        throw(ArgumentError("The number of src and destination knots must match."))
    end

    if length(src_knots) < 2 || length(dst_knots) < 2
        throw(ArgumentError("You need to specify at least two knots (the start and end of an interval)."))
    end

    # Promote the src_knots and dst_knots to a common type to address the case
    # where someone might use a mixture of integers and floating point numbers
    # when specifying the knots. 
    T₁ = promote_type(typeof.(src_knots)...)
    T₂ = promote_type(typeof.(dst_knots)...)
    T = promote_type(T₁, T₂)

    if saturate 
        intervals = convert_knots_to_intervals_and_saturate(T.(src_knots), T.(dst_knots))
    else
        intervals = convert_knots_to_intervals(T.(src_knots), T.(dst_knots))
    end   

    return PiecewiseLinearStretching(tuple(intervals...))
end

function convert_knots_to_intervals(src_knots::NTuple{N, <:Union{Real,AbstractGray}}, dst_knots::NTuple{N, <:Union{Real,AbstractGray}}) where {N}
    intervals = Vector{Pair{<:AbstractInterval, <: AbstractInterval}}(undef, N-1)
    for n = 1:(N-1)
        iₛ = src_knots[n] 
        jₛ = dst_knots[n] 
        iₑ = src_knots[n+1]
        jₑ = dst_knots[n+1]
        if n == 1
            interval = ClosedInterval(iₛ, iₑ) => ClosedInterval(jₛ, jₑ)
            intervals[n] = interval
        else
            interval = Interval{:open, :closed}(iₛ, iₑ) => Interval{:open, :closed}(jₛ, jₑ)
            intervals[n] = interval
        end
    end
    return intervals
end

function convert_knots_to_intervals_and_saturate(src_knots::NTuple{N, <:Union{Real,AbstractGray}}, dst_knots::NTuple{N, <:Union{Real,AbstractGray}}) where {N}
    intervals = convert_knots_to_intervals(src_knots, dst_knots)
    lb = gray(typemin(eltype(src_knots)))
    ub = gray(typemax(eltype(src_knots)))
    lb = lb == -Inf ? lb = nextfloat(-Inf) : lb
    ub = ub == Inf ? ub = prevfloat(Inf) : ub

    # Values that fall outside the first and last source knot get clamped to the first and
    # last destination knot, respectively. 
    interval₁ = Interval{:closed, :open}(lb, first(src_knots)) => ClosedInterval(first(dst_knots), first(dst_knots))
    interval₂ = Interval{:open, :closed}(last(src_knots), ub) => ClosedInterval(last(dst_knots), last(dst_knots))
    push!(intervals, interval₁)
    push!(intervals, interval₂)
    return intervals
end

function (f::PiecewiseLinearStretching)(out::GenericGrayImage, img::GenericGrayImage) 
    # Verify that the specified destination intervals fall inside the permissible range
    # of the output type. 
    T = eltype(out)
    for (_, dst) in f.intervals
        (a,b) = endpoints(dst)
        p = typemin(T)
        r = typemax(T)
        if  !(p <= a <= r) || !(p <= b <= r)
            throw(DomainError("The specified interval [$a, $b] falls outside the permissible range [$p, $r] associated with $T"))
        end
    end
    # A linear mapping of a value x in the interval [A, B] to the interval [a,b] is given by
    # the formula f(x) = (x-A) * ((b-a)/(B-A)) + a. The operation r * x - o is equivalent to
    # (x-A) * ((b-a)/(B-A)) + a. Precalculate the r and o so that later on an inner loop
    # contains only multiplication and addition to get better performance.
    N = length(f.intervals)
    rs = zeros(Float64, N)
    os = zeros(Float64, N)
    for (n, (src, dst)) in enumerate(f.intervals) 
        (A,B) = endpoints(src)
        (a,b) = endpoints(dst)   
        rs[n] = (b - a) / (B - A)
        os[n] = (A*b - B*a) / (B - A)
    end

    # Determine which interval the source pixel fall into, and apply the accompanying linear
    # transformation taking into account NaN values. 
    @inbounds for p in eachindex(img)
        val = img[p]
        if isnan(val)
            out[p] = val
        else
            is_inside_src_intervals = false
            for (n, (src, _)) in enumerate(f.intervals)                            
                if val ∈ src     
                    newval = rs[n] * val - os[n]                                                                          
                    out[p] = T <: Integer ? round(Int, newval) : newval  
                    is_inside_src_intervals = true  
                    break                                                                                   
                end                                                                                                    
            end 
            if !is_inside_src_intervals
                out[p] = val
            end
        end
    end
   return out                                                                                                                                                                                                             
end 

function (f::PiecewiseLinearStretching)(out::AbstractArray{<:Color3}, img::AbstractArray{<:Color3})
    T = eltype(out)
    yiq = convert.(YIQ, img)
    yiq_view = channelview(yiq)
    adjust_histogram!(view(yiq_view,1,:,:), f)
    out .= convert.(T, yiq)
end                                                                                                                                                                                               
