# Naive arbitrary resampler
# Uses a polyphase interpolator, followed by lineary interpolation
function naivefilt{T}( h::Vector{T}, x::Vector{T}, resamplerate::FloatingPoint, numfilters::Integer = 32 )

    xLen  = length( x )
    xx    = filt( h, x, numfilters//1 )
    xxLen = length(xx)
    yLen  = iceil( xLen * resamplerate )
    y     = similar( x, yLen )

    yIdx         = 1
    xxIdxVirtual = 1.0
    xxIdxLower   = 1
    xxIdxUpper   = 1
    Δ            = 0.0

    while xxIdxUpper <= xxLen

        println( "yIdx = $yIdx, xIdx = $(int(xxIdxVirtual/numfilters) + 1), xxIdxLower = $xxIdxLower of $xxLen, xxIdxUpper = $xxIdxUpper, Δ = $Δ" )

        y[yIdx] = xx[xxIdxLower] + Δ*( xx[xxIdxUpper] - xx[xxIdxLower] )

        yIdx         = yIdx + 1
        xxIdxVirtual = numfilters * ( yIdx - 1 ) / resamplerate + 1
        xxIdxLower   = ifloor( xxIdxVirtual )
        xxIdxUpper   = iceil( xxIdxVirtual )
        Δ            = xxIdxVirtual - xxIdxLower
    end

    resize!( y, yIdx-1 )

    return y
end
