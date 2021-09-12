using Base: depwarn

function adjust_histogram(img::Union{GenericGrayImage, AbstractArray{<:Color3}},
                           f::Union{LinearStretching, GammaCorrection},
                           args...; kwargs...)
    algo = typeof(f)
    depwarn("adjust_histogram(img, $algo) is deprecated, use adjust_intensity(img, $algo) instead", :adjust_histogram)                       
    return adjust_intensity(img, f, args...;kwargs...)
end

function adjust_histogram(type::Type{T},
                          img,
                          f::Union{LinearStretching, GammaCorrection},
                          args...; kwargs...)  where T
    algo = typeof(f)
    depwarn("adjust_histogram(::Type{T}, img, $algo) is deprecated, use adjust_intensity(::Type{T}, img, $algo) instead", :adjust_histogram)                       
    return adjust_intensity(type, img, f, args...;kwargs...)                    
end

function adjust_histogram(img::AbstractArray{T},
                          f::Union{LinearStretching, GammaCorrection},
                          args...; kwargs...) where T <: Colorant 
    algo = typeof(f)
    depwarn("adjust_histogram!(img, $algo) is deprecated, use adjust_intensity(img, $algo) instead", :adjust_histogram)                       
    return adjust_intensity(img, f, args...; kwargs...)
end

function adjust_histogram(type::Type{T},
                          img_sequence::Vector{<:AbstractArray},
                          f::Union{LinearStretching, GammaCorrection},
                          args...; kwargs...) where T
    algo = typeof(f)
    depwarn("adjust_histogram!(::Type{T}, img_sequence, $algo) is deprecated, use adjust_intensity(::Type{T}, img_sequence, $algo) instead", :adjust_histogram)       
    return adjust_histogram!(type, img_sequence, f, args...; kwargs...)               
end

function adjust_histogram!(img::Union{GenericGrayImage, AbstractArray{<:Color3}},
                           f::Union{LinearStretching, GammaCorrection},
                           args...; kwargs...)
    algo = typeof(f)
    depwarn("adjust_histogram!(img, $algo) is deprecated, use adjust_intensity!(img, $algo) instead", :adjust_histogram!)                       
    return adjust_intensity!(img, f, args...; kwargs...)                          
end

function adjust_histogram!(out_sequence::Vector{T},
                           img_sequence,
                           f::Union{LinearStretching, GammaCorrection},
                           args...; kwargs...) where T <: Union{GenericGrayImage, AbstractArray{<:Color3}} 
    algo = typeof(f)
    depwarn("adjust_histogram!(out_sequence, img_sequence, $algo) is deprecated, use adjust_intensity!(out_sequence, img_sequence, $algo) instead", :adjust_histogram!)
    return adjust_intensity(out_sequence, img_sequence, f, args...; kwargs...)
end