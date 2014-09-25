import Multirate
import Multirate: NaiveResamplers

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

numFilters   = 10
x            = [1.0:101]
resampleRate = 0.27

cutoffFreq      = 0.45
transitionWidth = 0.05
(hLen, β)       = Multirate.kaiserlength( transitionWidth, samplerate = numFilters )
hLen            = iceil(  hLen/numFilters  ) .* numFilters
h               = Multirate.firdes( hLen, cutoffFreq, DSP.kaiser, samplerate = numFilters, beta = β ) .* numFilters
h               = convert( Vector{Float32}, h)

@printf( "\n\tStateless arbitrary resampling\n\t\t" )
self = Multirate.FIRFilter( h, resampleRate, numFilters )
@time yStateless = Multirate.filt( self, x )

display( yStateless )
