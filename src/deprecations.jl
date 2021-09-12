using Base: depwarn

function adjust_histogram(img::Union{GenericGrayImage, AbstractArray{<:Color3}},
                           f::LinearStretching,
                           args...; kwargs...)
    
    depwarn("adjust_histogram(img, LinearStretching()) is deprecated, use adjust_intensity(img, LinearStretching()) instead", :adjust_histogram)                       
    return adjust_intensity(img, f, args...;kwargs...)
end

function adjust_histogram(type::Type{T},
                          img,
                          f::LinearStretching,
                          args...; kwargs...)  where T
    
    depwarn("adjust_histogram(::Type{T}, img, LinearStretching()) is deprecated, use adjust_intensity(::Type{T}, img, LinearStretching()) instead", :adjust_histogram)                       
    return adjust_intensity(type, img, f, args...;kwargs...)                    
end

function adjust_histogram(img::AbstractArray{T},
                          f::LinearStretching,
                          args...; kwargs...) where T <: Colorant 
    depwarn("adjust_histogram!(img, LinearStretching()) is deprecated, use adjust_intensity(img, LinearStretching()) instead", :adjust_histogram)                       
    return adjust_intensity(img, f, args...; kwargs...)
end

function adjust_histogram(type::Type{T},
                          img_sequence::Vector{<:AbstractArray},
                          f::LinearStretching,
                          args...; kwargs...) where T
    
    depwarn("adjust_histogram!(::Type{T}, img_sequence, LinearStretching()) is deprecated, use adjust_intensity(::Type{T}, img_sequence, LinearStretching()) instead", :adjust_histogram)       
    return adjust_histogram!(type, img_sequence, f, args...; kwargs...)               
end

function adjust_histogram!(img::Union{GenericGrayImage, AbstractArray{<:Color3}},
                           f::LinearStretching,
                           args...; kwargs...)

    depwarn("adjust_histogram!(img, LinearStretching()) is deprecated, use adjust_intensity!(img, LinearStretching()) instead", :adjust_histogram!)                       
    return adjust_intensity!(img, f, args...; kwargs...)                          
end

function adjust_histogram!(out_sequence::Vector{T},
                           img_sequence,
                           f::LinearStretching,
                           args...; kwargs...) where T <: Union{GenericGrayImage, AbstractArray{<:Color3}} 
    
    depwarn("adjust_histogram!(out_sequence, img_sequence, LinearStretching()) is deprecated, use adjust_intensity!(out_sequence, img_sequence, LinearStretching()) instead", :adjust_histogram!)
    return adjust_intensity(out_sequence, img_sequence, f, args...; kwargs...)
end