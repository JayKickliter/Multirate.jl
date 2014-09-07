# Naive arbitrary resampler
# Uses a polyphase interpolator, followed by lineary interpolation
function naivefilt{T}( h::Vector{T}, x::Vector{T}, resamplerate::FloatingPoint, numfilters::Integer = 32 )
    xLen = length( x )
    yLen = iceil( xLen * resamplerate )
    y    = similar( x, yLen )
    yIdx = 1
    xx   = filt( h, x, numfilters//1 )

    for yIdx in 1:yLen
        xxIdxVirtual = numfilters * ( yIdx - 1 ) / resamplerate + 1
        xxIdxLower   = ifloor( xxIdxVirtual )
        xxIdxUpper   = iceil( xxIdxVirtual )
        Δ            = xxIdxVirtual - xxIdxLower
        y[yIdx]      = xx[xxIdxLower] + Δ*( xx[xxIdxUpper] - xx[xxIdxLower] )
    end

    return y
end