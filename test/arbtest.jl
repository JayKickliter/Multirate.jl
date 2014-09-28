import Multirate
import Multirate: NaiveResamplers
using Base.Test

function Base.isapprox( x1::Vector, x2::Vector )
    Nx1 = length( x1 )
    Nx2 = length( x2 )

    if Nx1 != Nx2
        @printf("x1 & x2 are different lengths vectors\n")
        return false
    end

    for i = 1:Nx1
        if !isapprox( x1[i], x2[i] )
            @printf( "Something went wrong at index %d\n", i )
            return false
        end
    end

    return true
end

Tx           = Float32
Th           = Float32
numFilters   = 32
x            = Tx[1.0:101]
resampleRate = 33/32

cutoffFreq      = 0.45
transitionWidth = 0.05
(hLen, β)       = Multirate.kaiserlength( transitionWidth, samplerate = numFilters )
hLen            = iceil(  hLen/numFilters  ) .* numFilters
h               = Multirate.firdes( hLen, cutoffFreq, DSP.kaiser, samplerate = numFilters, beta = β ) .* numFilters
h               = convert( Vector{Th}, h)

@time yNaive = NaiveResamplers.naivefilt( h, x, resampleRate, numFilters )
@time yArb   = Multirate.filt( h, x, resampleRate )

commonLen = min( length(yNaive), length(yArb) )
isapprox( yNaive[1:commonLen], yArb[1:commonLen] ) ||  display( [ [1:commonLen] yNaive[1:commonLen] yArb[1:commonLen] yNaive[1:commonLen].-yArb[1:commonLen]] )
