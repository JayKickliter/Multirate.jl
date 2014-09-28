using Winston
using Polynomials
using Multirate




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

function polyfit( y::AbstractVector, polyorder::Integer )
  A = [ x^p for x in 1:length(y), p = 0:polyorder ]
  Poly(A \ y)
end

# Convert a polyphase filterbank into a polynomial filterbank
function pfb2pnfb{T}( pfb::Matrix{T}, polyorder = 5 )
    (tapsPer𝜙, N𝜙) = size( pfb )
    result         = Array( Poly{T}, tapsPer𝜙 )
    sizehint( result, tapsPer𝜙 )
    
    for i in 1:tapsPer𝜙
        row = vec( pfb[i,:] )
        result[i] = polyfit( row, polyorder )
    end 

    return result
end

function pnfb2pfb{T}( pnfb::Vector{Poly{T}}, N𝜙::Integer )
    tapsPer𝜙 = length(pnfb)
    pfb      = Array( T, tapsPer𝜙, N𝜙 )
    
    for i in 1:N𝜙, j in 1:tapsPer𝜙
        pfb[j,i] = polyval( pnfb[j], i )
    end
    
    return pfb
end

function test()
    numFilters      = 10
    resampleRate    = 0.98
    cutoffFreq      = 0.4
    transitionWidth = 0.1
    (hLen, β)       = kaiserlength( transitionWidth, samplerate = numFilters )
    hLen            = iceil(  hLen/numFilters  ) .* numFilters
    h               = Multirate.firdes( hLen, cutoffFreq, DSP.kaiser, samplerate = numFilters, beta = β ) 
    pfb             = taps2pfb( h, numFilters )
    (tapsPer𝜙, N𝜙)  = size( pfb )
    polyorder       = 6
    
    pnfb = pfb2pnfb( pfb, polyorder )
    
    return pnfb, pfb
end

(pnfb, pfb) = test()
pfbr        = pnfb2pfb( pnfb, 10 )
absdpfb     = abs(pfb.-pfbr)
x           = [1:]
# display(pnfb)
# display(pfb)
# display(pfbr)
display( absdpfb )
display( (maximum(absdpfb), minimum(absdpfb)) )
plot( )
