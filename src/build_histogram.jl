function partition_interval(nbins::Integer, minval::Real, maxval::Real)
    edges = (range(0,step=(maxval - minval) / nbins,length=nbins)) .+ minval
end

function partition_interval(nbins::Integer, minval::AbstractGray, maxval::AbstractGray)
    partition_interval(nbins, minval.val, maxval.val)
end

"""
```
edges, count = build_histogram(img, nbins)
edges, count = build_histogram(img, nbins; minval, maxval)
edges, count = build_histogram(img, edges)
```

Generates a histogram for the image over `nbins` spread between `[minval, maxval]`.
Color images are automatically converted to grayscale.

# Output
Returns `edges` which is a `AbstractRange` type that specifies how the  interval
`[minval, maxval]` is divided into bins, and an array `count` which records the
concomitant bin frequencies. In particular, `count` has the following
properties:

* `count[0]` is the number satisfying `x < edges[1]`
* `count[i]` is the number of values `x` that satisfy `edges[i] <= x < edges[i+1]`
* `count[end]` is the number satisfying `x >= edges[end]`.
* `length(count) == length(edges)+1`.

# Details

One can consider a histogram as a piecewise-constant model of a probability
density function ``f`` [1]. Suppose that ``f`` has support on some interval ``I =
[a,b]``.  Let ``m`` be an integer and ``a = a_1 < a_2 < \\ldots < a_m < a_{m+1} =
b`` a sequence of real numbers. Construct a sequence of intervals

```math
I_1 = [a_1,a_2], I_2 = (a_2, a_3], \\ldots, I_{m} = (a_m,a_{m+1}]
```

which partition ``I`` into subsets ``I_j`` ``(j = 1, \\ldots, m)`` on which
``f`` is constant. These subsets satisfy ``I_i \\cap I_j = \\emptyset, \\forall
i \\neq j``, and are commonly referred to as *bins*. Together they encompass the
entire range of data values such that ``\\sum_j |I_j | = | I |``. Each bin has
width ``w_j = |I_j| = a_{j+1} - a_j`` and height ``h_j`` which is the constant
probability density over the region of the bin. Integrating the constant
probability density over the width of the bin ``w_j`` yields a probability mass
of ``\\pi_j = h_j w_j`` for the bin.

For a sample ``x_1, x_2, \\ldots, x_N``, let


```math
n_j = \\sum_{n = 1}^{N}\\mathbf{1}_{(I_j)}(x_n),
\\quad \\text{where} \\quad
\\mathbf{1}_{(I_j)}(x) =
\\begin{cases}
 1 & \\text{if} \\; x \\in I_j,\\\\
 0 & \\text{otherwise},
\\end{cases},
```
represent the number of samples falling into the interval ``I_j``. An estimate
for the probability mass of the ``j``th bin is given by the relative frequency
``\\hat{\\pi} = \\frac{n_j}{N}``, and the histogram estimator of the probability
density function is defined as
```math
\\begin{aligned}
\\hat{f}_n(x)  & = \\sum_{j = 1}^{m}\\frac{n_j}{Nw_j} \\mathbf{1}_{(I_j)}(x) \\\\
& = \\sum_{j = 1}^{m}\\frac{\\hat{\\pi}_j}{w_j} \\mathbf{1}_{(I_j)}(x) \\\\
& = \\sum_{j = 1}^{m}\\hat{h}_j \\mathbf{1}_{(I_j)}(x).
\\end{aligned}
```

The function ``\\hat{f}_n(x)`` is a genuine density estimator because ``\\hat{f}_n(x)  \\ge 0`` and
```math
\\begin{aligned}
\\int_{-\\infty}^{\\infty}\\hat{f}_n(x) \\operatorname{d}x & = \\sum_{j=1}^{m} \\frac{n_j}{Nw_j} w_j \\\\
& = 1.
\\end{aligned}
```

# Options
Various options for the parameters of this function are described in more detail
below.

## Choices for `nbins`
You can specify the number of discrete bins for the histogram. When specifying
the number of bins consider the maximum number of graylevels that your image
type supports. For example, with an image of type `N0f8` there is a maximum
of 256 possible graylevels. Hence, if you request more than 256 bins for
that type of image you should expect to obtain zero counts for numerous bins.


## Choices for `minval`
You have the option to specify the lower bound of the interval over which the
histogram will be computed.  If `minval` is not specified then the minimum
value present in the image is taken as the lower bound.

## Choices for `maxval`
You have the option to specify the upper bound of the interval over which the
histogram will be computed.  If `maxval` is not specified then the maximum
value present in the image is taken as the upper bound.

## Choices for `edges`
If you do not designate the number of bins, nor the lower or upper bound of the
interval, then you have the option to directly stipulate how the intervals will
be divided by specifying a `AbstractRange` type.

# Example

Compute the histogram of a grayscale image.
```julia

using TestImages, FileIO, ImageView

img =  testimage("mandril_gray");
edges, counts  = build_histogram(img, 256, minval = 0, maxval = 1)
```

Given a color image, compute the histogram of the red channel.
```julia
img = testimage("mandrill")
r = red.(img)
edges, counts  = build_histogram(r, 256, minval = 0, maxval = 1)
```

# References
[1] E. Herrholz, "Parsimonious Histograms," Ph.D. dissertation, Inst. of Math. and Comp. Sci., University of Greifswald, Greifswald, Germany, 2011.
"""
function build_histogram(img::AbstractArray, nbins::Integer = 256; minval::Union{Real,AbstractGray}, maxval::Union{Real,AbstractGray})
    edges = partition_interval(nbins, minval, maxval)
    build_histogram(img, edges)
end

# Warning is due to: https://github.com/JuliaLang/julia/issues/9498
# Might resolve this based on this idea: https://discourse.julialang.org/t/fun-with-kwargs-methods/15711/8
function build_histogram(img::AbstractArray, nbins::Integer = 256)
    build_histogram(img, nbins, minval =  minfinite(img), maxval = maxfinite(img))
end

function build_histogram(img::AbstractArray{T}, edges::AbstractRange) where {T<:Color3}
    build_histogram(Gray.(img),edges)
end

function build_histogram(img::AbstractArray, edges::AbstractRange)
    Base.has_offset_axes(edges) && throw(ArgumentError("edges must be indexed starting with 1"))
    lb = first(axes(edges,1))-1
    ub = last(axes(edges,1))
    first_edge = first(edges)
    step_size = step(edges)
    counts = fill(0, lb:ub)
    for val in img
         if isnan(val)
             continue
         else
            if val >= edges[end]
                counts[end] += 1
                continue
            elseif val < first(edges)
                counts[lb] += 1
            else
                index = Int(Base.div(val-first_edge,step_size)) + 1
                counts[index] += 1
            end
        end
    end
    edges, counts
end
