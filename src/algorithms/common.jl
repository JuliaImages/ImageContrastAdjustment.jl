function transform_density!(img::GenericGrayImage, edges::AbstractRange, cdf::AbstractVector, minval::Union{Real,AbstractGray}, maxval::Union{Real,AbstractGray})
    first_edge, last_edge = first(edges), last(edges)
    first_cdf, last_cdf = first(cdf), last(cdf)
    inv_step_size = 1/step(edges)
    scale = (maxval - minval) / (cdf[end] - first(cdf))
    function transform(val)
        val = gray(val)
        if val >= last_edge
            newval = last_cdf
        elseif val < first_edge
            newval = first_cdf
        else
            index = floor(Int, (val-first_edge)*inv_step_size) + 1
            @inbounds newval = cdf[index]
        end
        # Scale the new intensity value to so that it lies in the range [minval, maxval].
        newval = minval + (newval - first_cdf) * scale
    end
    if eltype(img) <: Integer
        map!(val->ceil(transform(val)), img, img)
    else
        map!(transform, img, img)
    end
end

function construct_pdfs(img::GenericGrayImage, targetimg::AbstractArray, edges::AbstractRange)
    _, histogram = build_histogram(img, edges)
    _, target_histogram = build_histogram(targetimg, edges)
    return edges, histogram / sum(histogram), target_histogram / sum(target_histogram)
end

function construct_pdfs(img::GenericGrayImage, targetimg::AbstractArray, nbins::Integer = 256)
    if eltype(img) <: AbstractGray
        imin, imax = 0, 1
    else
        imin, imax = min(minfinite(img), minfinite(targetimg)), max(maxfinite(img), maxfinite(targetimg))
    end
    edges, histogram = build_histogram(img, nbins, minval = imin, maxval = imax)
    _, target_histogram = build_histogram(targetimg, edges)
    return edges, histogram / sum(histogram), target_histogram / sum(target_histogram)
end

function lookup_icdf(cdf::AbstractArray, targetcdf::AbstractArray)
    lookup_table = zeros(Int, length(cdf))
    i = 1
    for j = 1:length(cdf)
        p = cdf[j]
        while i < length(targetcdf) && targetcdf[i+1] <= p
            i += 1
        end
        lookup_table[j] = i
    end
    lookup_table
end

function cdf2pdf!(pdf::AbstractArray, cdf::AbstractArray)
    pdf[1] = cdf[1]
    for i = length(cdf)-1:-1:2
        pdf[i] = cdf[i] - cdf[i-1]
    end
end
