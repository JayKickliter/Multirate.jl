module NaiveResamplers
export naivefilt

# Naive rational resampler
function naivefilt{T}( h::Vector{T}, x::Vector{T}, resamplerate::Rational = 1//1 )

    upfactor     = num( resamplerate )
    downfactor   = den( resamplerate )
    xLen         = length( x )
    xZeroStuffed = zeros( T, length( x ) * upfactor )
    yLen         = int( ceil( length(x) * resamplerate ) )
    y            = Array( T, yLen )

    for n in 0:length(x)-1
        xZeroStuffed[ n*upfactor+1 ] = x[ n+1 ]
    end

    y = firfilt( h, xZeroStuffed )
    y = [ y[n] for n = 1:downfactor:length( y ) ]
end




# Naive arbitrary resampler
# Uses a polyphase interpolator, followed by lineary interpolation
function naivefilt{T}( h::Vector{T}, x::Vector{T}, resamplerate::FloatingPoint, numfilters::Integer = 32 )

    xLen          = length( x )
    xInterpolated = naivefilt( h, x, numfilters//1 )
    xxLen         = length(xInterpolated)
    yLen          = iceil( xLen * resamplerate )
    y             = similar( x, yLen )
    yIdx          = 1
    xxIdxVirtual  = 1.0
    xxIdxLower    = 1
    xxIdxUpper    = 1
    Δ             = 0.0

    while xxIdxUpper <= xxLen

        # println( "yIdx = $yIdx, xIdx = $(int(xxIdxVirtual/numfilters) + 1), xxIdxLower = $xxIdxLower of $xxLen, xxIdxUpper = $xxIdxUpper, Δ = $Δ" )

        y[yIdx] = xInterpolated[xxIdxLower] + Δ*( xInterpolated[xxIdxUpper] - xInterpolated[xxIdxLower] )

        yIdx         = yIdx + 1
        xxIdxVirtual = numfilters * ( yIdx - 1 ) / resamplerate + 1
        xxIdxLower   = ifloor( xxIdxVirtual )
        xxIdxUpper   = iceil( xxIdxVirtual )
        Δ            = xxIdxVirtual - xxIdxLower
    end

    resize!( y, yIdx-1 )

    return y
end

end # module