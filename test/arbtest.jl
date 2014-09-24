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
# x            = [1.0:101]
x = ones(101)
resampleRate = 0.27

cutoffFreq      = 0.45
transitionWidth = 0.05
(hLen, β)       = Multirate.kaiserlength( transitionWidth, samplerate = numFilters )
hLen            = iceil(  hLen/numFilters  ) .* numFilters
h               = Multirate.firdes( hLen, cutoffFreq, DSP.kaiser, samplerate = numFilters, beta = β ) .* numFilters
# h               = ones(numFilters)./numFilters

@printf( "\n\tNaive arbitrary resampling\n\t\t" )
@time naiveResult = NaiveResamplers.naivefilt( h, x, resampleRate, numFilters )

@printf( "\n\tStateless arbitrary resampling\n\t\t" )
@time statelessResult = Multirate.filt( h, x, resampleRate, numFilters )

@printf( "\n\tPiecewise arbitrary resampling\n\t\t" )
self           = Multirate.FIRFilter( h, resampleRate, numFilters )
piecwiseResult = eltype(x)[]
sizehint( piecwiseResult, iceil( length(x)*resampleRate ) )
@time for i in 1:length( x )
    thisY = filt( self, x[i:i] )
    append!( piecwiseResult, thisY )
end

commonLen = min( length(naiveResult), length(statelessResult), length(piecwiseResult) )

resize!( naiveResult, commonLen )
resize!( statelessResult, commonLen )
resize!( piecwiseResult, commonLen )

display( [  [1:commonLen] naiveResult statelessResult  piecwiseResult abs(naiveResult.-statelessResult) abs(naiveResult.-piecwiseResult) ] )
display( [ [1:commonLen-1] diff(naiveResult) diff(statelessResult) diff(piecwiseResult) ] )