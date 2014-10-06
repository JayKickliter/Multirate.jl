# Interpolator FIR kernel
type Channelizer{T} <: FIRKernel
    pfb::PFB{T}
    Nchannels::Int
    tapsPer𝜙::Int
    history::Vector
end

function Channelizer( h::Vector, Nchannels::Integer )
    pfb       = taps2pfb( h, Nchannels )
    Nchannels = size( pfb )[2]
    tapsPer𝜙  = size( pfb )[1]
    Channelizer( pfb, Nchannels, tapsPer𝜙, [] )
end

function Channelizer( Nchannels::Integer, tapsPer𝜙 = 10 )
    hLen = tapsPer𝜙 * Nchannels
    h    = firdes( hLen, 0.45/Nchannels, kaiser ) .* Nchannels
    Channelizer( h, Nchannels )
end



function filt!{Tb,Th,Tx}( buffer::Matrix{Tb}, kernel::Channelizer{Th}, x::Vector{Tx} )
    xLen                = length( x )
    (bufLen,bufWidth)   = size( buffer )
    fftBuffer           = Array( Tb, kernel.Nchannels )

    @assert xLen % kernel.Nchannels == 0
    @assert bufLen * bufWidth       == xLen
    
    if kernel.history == []
        kernel.history = Array(Vector, kernel.Nchannels)
        for channel in 1:kernel.Nchannels
            kernel.history[channel] = zeros(Tx, kernel.tapsPer𝜙-1)
        end
    end
    
    # TODO: splitup x similar to polyize

    𝜙Idx   = kernel.Nchannels
    xIdx   = 1
    rowIdx = 1
    while xIdx <= xLen

        if xIdx < kernel.tapsPer𝜙
            history = kernel.history[𝜙Idx]
            fftBuffer[𝜙Idx] = unsafedot( kernel.pfb, 𝜙Idx, history, x, xIdx )
        else
            fftBuffer[𝜙Idx] = unsafedot( kernel.pfb, 𝜙Idx, x, xIdx )
        end
        
        xIdx += 1
        𝜙Idx -= 1
        
        if 𝜙Idx == 0
            buffer[rowIdx,:] = flipud(fft(fftBuffer))
            𝜙Idx             = kernel.Nchannels
            rowIdx          += 1
        end        
    end
    
    return buffer
end

function filt{Th,Tx}( kernel::Channelizer{Th}, x::Vector{Tx} )
    xLen   = length( x )

    @assert xLen % kernel.Nchannels == 0

    buffer = Array( promote_type(Th,Tx), int(xLen/kernel.Nchannels), kernel.Nchannels )
    filt!( buffer, kernel, x )
    return buffer
end
