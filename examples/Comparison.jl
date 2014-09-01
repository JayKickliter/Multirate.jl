import DSP: firfilt
import Multirate: filt, firdes

sampleRate    = 48000
interpolation = 147
decimation    = 160
ratio         = interpolation  // decimation
numTaps       = 24*interpolation
x             = rand( 1_000_000 )
h             = firdes( numTaps, 0.5/interpolation, Multirate.kaiser, beta = 7.8562  )

function bruteresample{T}( h::Vector{T}, x::Vector{T}, ratio::Rational{Int} )
    upfactor   = num( ratio )
    downfactor = den( ratio )
    xStuffed   = zeros( T, length(x) * upfactor )

    for n = 0:length(x)-1;
        xStuffed[ n*upfactor+1 ] = x[ n+1 ]
    end

    yInterpolated = firfilt( h, xStuffed )
    y = [ yInterpolated[n] for n = 1:downfactor:length( yInterpolated ) ]
end

@time y = bruteresample( h, x, 147//160 );
@time y = filt( h, x, ratio );
