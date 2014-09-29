using Polynomials

typealias PNFB{T} Vector{Poly{T}}




function polyfit( y::AbstractVector, polyorder::Integer )
  A = [ x^p for x in 1:length(y), p = 0:polyorder ]
  Poly(A \ y)
end




# Convert a polyphase filterbank into a polynomial filterbank
function pfb2pnfb{T}( pfb::PFB{T}, polyorder )
    (tapsPer𝜙, N𝜙) = size( pfb )
    result         = Array( Poly{T}, tapsPer𝜙 )

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

type FIRFarrow{T} <: FIRKernel
    rate::Float64
    pfb::PFB{T}
    pnfb::PNFB{T}
    currentTaps::Vector{T}
    N𝜙::Int
    tapsPer𝜙::Int
    𝜙Idx::Float64
    Δ::Float64
    inputDeficit::Int
    xIdx::Int
end




function FIRFarrow{Th}( h::Vector{Th}, rate::Real, N𝜙::Integer, polyorder::Integer )
    pfb          = flipud(taps2pfb( h,  N𝜙 ))
    pnfb         = pfb2pnfb( pfb, polyorder )
    tapsPer𝜙     = size( pfb )[1]
    𝜙Idx         = 1.0
    Δ            = N𝜙/rate
    inputDeficit = 1
    xIdx         = 1
    currentTaps  = Th[ polyval( pnfb[tapIdx], 𝜙Idx ) for tapIdx in 1:tapsPer𝜙 ]
    FIRFarrow( rate, pfb, pnfb, currentTaps, N𝜙, tapsPer𝜙, 𝜙Idx, Δ, inputDeficit, xIdx )
end




function update!( kernel::FIRFarrow )
    kernel.𝜙Idx += kernel.Δ

    if kernel.𝜙Idx > kernel.N𝜙
        kernel.xIdx += ifloor( (kernel.𝜙Idx-1) / kernel.N𝜙 )
        kernel.𝜙Idx  = mod( (kernel.𝜙Idx-1), kernel.N𝜙 ) + 1
    end

    for tapIdx in 1:kernel.tapsPer𝜙
        @inbounds kernel.currentTaps[tapIdx] = polyval( kernel.pnfb[tapIdx], kernel.𝜙Idx )
    end
end




function filt{Th,Tx}( self::FIRFilter{FIRFarrow{Th}}, x::Vector{Tx} )
    kernel              = self.kernel
    xLen                = length( x )
    bufLen              = iceil( xLen * kernel.rate ) + 1
    buffer              = zeros(promote_type(Th,Tx), bufLen)
    bufIdx              = 1
    history::Vector{Tx} = self.history
    db_vec_phi          = Array(Float64, bufLen)
    db_vec_xidx         = Array(Int, bufLen)

    # Do we have enough input samples to produce one or more output samples?
    if xLen < kernel.inputDeficit
        self.history = shiftin!( history, x )
        kernel.inputDeficit -= xLen
        return buffer[1:bufIdx-1]
    end

    # Skip over input samples that are not needed to produce output results.
    kernel.xIdx = kernel.inputDeficit

    while kernel.xIdx <= xLen
        db_vec_xidx[bufIdx] = kernel.xIdx
        db_vec_phi[bufIdx]  = kernel.𝜙Idx
        if kernel.xIdx < kernel.tapsPer𝜙
            y = unsafedot( kernel.currentTaps, history, x, kernel.xIdx )
        else
            y = unsafedot( kernel.currentTaps, x, kernel.xIdx )
        end
        buffer[bufIdx] = y
        bufIdx        += 1
        update!( kernel )
    end

    # Did we overestimate needed buffer size?
    # TODO: Get rid of this by correctly calculating output size.
    bufLen == bufIdx - 1 || resize!( buffer, bufIdx - 1)
    kernel.inputDeficit = kernel.xIdx - xLen

    self.history = shiftin!( history, x )

    resize!( db_vec_phi, length(buffer) )
    resize!( db_vec_xidx, length(buffer) )    

    return buffer, db_vec_xidx, db_vec_phi
end




function FIRFilter( h::Vector, rate::FloatingPoint, N𝜙::Integer, polyorder::Integer )
    rate > 0.0 || error( "rate must be greater than 0" )
    kernel     = FIRFarrow( h, rate, N𝜙, polyorder )
    historyLen = kernel.tapsPer𝜙 - 1
    history    = zeros( historyLen )
    FIRFilter( kernel, history, historyLen )
end




function filt( h::Vector, x::Vector, rate::FloatingPoint, N𝜙::Integer, polyorder::Integer )
    self = FIRFilter( h, rate, N𝜙, polyorder )
    filt( self, x )
end