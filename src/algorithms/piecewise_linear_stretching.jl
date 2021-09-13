"""
```
    PiecewiseLinearStretching <: AbstractHistogramAdjustmentAlgorithm
    PiecewiseLinearStretching(;  src_intervals = (0.0 => 0.1, 0.1 => 0.9, 0.9 => 1.0), 
                                 dst_intervals = (0.0 => 0.2, 0.2 => 1.0, 1.0 => 1.0))

    adjust_histogram([T,] img, f::PiecewiseLinearStretching)
    adjust_histogram!([out,] img, f::PiecewiseLinearStretching)
```

Returns an image for which the intensity range was partitioned into subintervals (specified
by `src_intervals`) and mapped linearly to a new set of subintervals
(specified by `dst_intervals`). 

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

## Choices for `src_intervals` 

The `src_intervals` must be specified as an `ntuple` of pairs. For example,
`src_intervals = (0.0 => 0.2, 0.2 => 0.8, 0.8 => 1.0)` partitions the unit
interval into three subintervals, [0.0, 0.2], (0.2, 0.8] and (0.8, 1.0]. 

To specify a single interval you must include a trailing comma, e.g.
`src_intervals = (0.0 => 1.0,)`. 

The start and end of a subinterval in `src_intervals` cannot be equal. For example,
`src_intervals = (0.0 => 0.5, 0.5 => 0.5, 0.5 => 1.0)` will result in an `ArgumentError`
because of the pair `0.5 => 0.5`.

## Choices for `dst_intervals` 

The `dst_intervals` must be specified as an `ntuple` of pairs and there needs to be a
corresponding pair for each pair in `src_intervals`.
    
Unlike `src_intervals`, the start and end of an interval in `src_intervals` can be equal.
For example, the combination `src_intervals = (0.0 => 0.5, 0.5 => 1.0)` and `dst_intervals =
(0.0 => 0.0 , 1.0 => 1.0)` will produce a binary image. 

Care must be taken to ensure that the subintervals specified in `dst_intervals` are
representable within the type of image on which the piecewise linear stretching will be
applied. For example, if you specify a negative output interval like `dst_intervals = (1.0
=> -1.0),` and try to apply the transformation on an image with type `Gray{N0f8}`, you will
receive a DomainError. 

## Saturating pixels

If the specified `src_intervals` don't span the entire range of intensities in an image,
then intensities that fall outside `src_intervals` are automatically set to the boundary
values of the `dst_intervals`. For example, given the combination `src_intervals = (0.1 =>
0.9, )` and `dst_intervals = (0.0 => 1.0)`, intensities below `0.1` are set to zero, and
intensities above `0.9` are set to one. 

# Example

```julia
using ImageContrastAdjustment, TestImages

img = testimage("mandril_gray")
# Stretches the contrast in `img` so that it spans the unit interval.
f = PiecewiseLinearStretching(src_intervals = (0.0 => 0.2, 0.2 => 0.8, 0.8 => 1.0),
                              dst_intervals = (0.0 => 0.4, 0.4 => 0.9, 0.9 => 1.0))
imgo = adjust_histogram(img, f)
```

# References
1. W. Burger and M. J. Burge. *Digital Image Processing*. Texts in Computer Science, 2016. [doi:10.1007/978-1-4471-6684-9](https://doi.org/10.1007/978-1-4471-6684-9)
"""
@with_kw struct PiecewiseLinearStretching{T} <: AbstractHistogramAdjustmentAlgorithm where T <: NTuple{N, Pair{<: Union{Real, AbstractGray}, <: Union{Real, AbstractGray}}} where {N}
    src_intervals::T = (0.0 => 0.1, 0.1 => 0.9, 0.9 => 1.0)
    dst_intervals::T = (0.0 => 0.2, 0.2 => 1.0, 1.0 => 1.0)
    function PiecewiseLinearStretching(src_intervals::T₁, dst_intervals::T₂) where {T₁ <: NTuple{N₁, Pair{<: Union{Real, AbstractGray}, <: Union{Real, AbstractGray}}} where {N₁}, T₂ <: NTuple{N₂, Pair{<: Union{Real, AbstractGray}, <: Union{Real, AbstractGray}}} where {N₂}}  
        length(src_intervals) == length(dst_intervals) ||  throw(ArgumentError("The dst_intervals must define an interval pair for each consecutive src_interval pair, i.e. length(src_intervals) == length(dst_intervals)."))
        for src_pair in src_intervals
            if first(src_pair) == last(src_pair) 
                throw(ArgumentError("The start and end of an interval in src_intervals cannot be equal."))
            end
        end
        T = promote_type(T₁, T₂)
        new{T}(src_intervals, dst_intervals)
    end    
end

function (f::PiecewiseLinearStretching)(out::GenericGrayImage, img::GenericGrayImage) 
    # Verify that the specified dst_intervals fall inside the permissible range
    # of the output type. 
    T = eltype(out)
    for (a,b) in f.dst_intervals
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
    N = length(f.src_intervals)
    rs = zeros(Float64, N)
    os = zeros(Float64, N)
    for (src, dst) in zip(pairs(f.src_intervals), pairs(f.dst_intervals)) 
        n, (A,B) = src
        _, (a,b) = dst   
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
            for (n, (A,B)) in pairs(f.src_intervals)                                 
                if (A <= val < B) || (A <= val <= B && n == N)       
                    newval = rs[n] * val - os[n]                                                                          
                    out[p] = T <: Integer ? round(Int, newval) : newval  
                    is_inside_src_intervals = true  
                    break                                                                                   
                end                                                                                                    
            end 
            if !is_inside_src_intervals
                # Saturate pixels that fall outside the specified src_intervals by setting
                # them equal to the start/end of the edge dst_intervals. 
                a = first(first(f.dst_intervals))
                b = last(last(f.dst_intervals))
                A = first(first(f.src_intervals))
                newval = val < A ? a : b
                out[p] = T <: Integer ? round(Int, newval) : newval  
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
                                                                                                                                                                                                   
