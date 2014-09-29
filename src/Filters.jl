#==============================================================================#
#       ___ _   _ ___  ____ ____      /    ____ ____ _  _ ____ ___ ____        #
#        |   \_/  |__] |___ [__      /     |    |  | |\ | [__   |  |__/        #
#        |    |   |    |___ ___]    /      |___ |__| | \| ___]  |  |  \ .      #
#==============================================================================#

typealias PFB{T} Matrix{T}

abstract Filter
abstract FIRKernel

# Single rate FIR kernel, just hold filter h
type FIRStandard <: FIRKernel
    h::Vector
    hLen::Int
    function FIRStandard( h::Vector )
        self      = new()
        self.h    = flipud( h )
        self.hLen = length( h )
        return self
    end
end


# Interpolator FIR kernel
type FIRInterpolator <: FIRKernel
    pfb::PFB
    interpolation::Int
    N𝜙::Int
    tapsPer𝜙::Int
    function FIRInterpolator( h::Vector, interpolation::Integer )
        self               = new()
        self.pfb           = flipud( taps2pfb( h, interpolation ) )
        self.tapsPer𝜙      = size( self.pfb )[1]
        self.N𝜙            = size( self.pfb )[2]
        self.interpolation = interpolation
        return self
    end
end


# Decimator FIR kernel
type FIRDecimator <: FIRKernel
    h::Vector
    hLen::Int
    decimation::Int
    inputDeficit::Int
    function FIRDecimator( h::Vector, decimation::Integer )
        self              = new()
        self.h            = flipud( h )
        self.hLen         = length( h )
        self.decimation   = decimation
        self.inputDeficit = 1
        return self
    end
end


# Rational resampler FIR kernel
type FIRRational  <: FIRKernel
    pfb::PFB
    ratio::Rational{Int}
    N𝜙::Int
    tapsPer𝜙::Int
    criticalYidx::Int
    𝜙Idx::Int
    inputDeficit::Int
    function FIRRational( h::Vector, ratio::Rational )
        self              = new()
        self.pfb          = flipud( taps2pfb( h, num(ratio) ))
        self.ratio        = ratio
        self.N𝜙           = size( self.pfb )[2]
        self.tapsPer𝜙     = size( self.pfb )[1]
        self.criticalYidx = ifloor( self.tapsPer𝜙 * ratio )
        self.𝜙Idx         = 1
        self.inputDeficit = 1
        return self
    end
end


# Arbitrary resampler FIR kernel
# This kernel is different from the others in that it has two polyphase filtlter banks.
# The the second filter bank, dpfb, is the derivative of pfb. The purpose of this is to
# allow us to compute two y values, yLower & yUpper, whitout having to advance the input
# index by 1. It makes the kernel simple by not having to store extra state in the case
# when where's at the last polphase branch and the last available input sample. By using
# a derivitive filter, we can always compute the output in that scenario.
# See section 7.6.1 in [1] for a better explanation.

type FIRArbitrary{T} <: FIRKernel
    rate::Float64
    pfb::PFB{T}
    dpfb::PFB{T}
    N𝜙::Int
    𝜙Stride::Int
    tapsPer𝜙::Int
    𝜙Idx::Int
    α::Float64
    δ::Float64
    inputDeficit::Int
    xIdx::Int
end

function FIRArbitrary( h::Vector, rate::Real, N𝜙::Integer )
    dh           = [ diff( h ), 0 ]
    pfb          = flipud(taps2pfb( h,  N𝜙 ))
    dpfb         = flipud(taps2pfb( dh, N𝜙 ))
    tapsPer𝜙     = size( pfb )[1]
    𝜙Idx         = 1
    α            = 0.0
    (δ, 𝜙Stride) = modf( N𝜙/rate )
    inputDeficit = 1
    xIdx         = 1
    FIRArbitrary( rate, pfb, dpfb, N𝜙, int( 𝜙Stride ), tapsPer𝜙, 𝜙Idx, α, δ, inputDeficit, xIdx )
end



# FIRFilter - the kernel does the heavy lifting
type FIRFilter{Tk<:FIRKernel} <: Filter
    kernel::Tk
    history::Vector
    historyLen::Int
end

function FIRFilter( h::Vector, resampleRatio::Rational = 1//1 )
    interpolation = num( resampleRatio )
    decimation    = den( resampleRatio )
    historyLen    = 0

    if resampleRatio == 1                                     # single-rate
        kernel     = FIRStandard( h )
        historyLen = kernel.hLen - 1
    elseif interpolation == 1                                 # decimate
        kernel     = FIRDecimator( h, decimation )
        historyLen = kernel.hLen - 1
    elseif decimation == 1                                    # interpolate
        kernel        = FIRInterpolator( h, interpolation )
        historyLen = kernel.tapsPer𝜙 - 1
    else                                                      # rational
        kernel     = FIRRational( h, resampleRatio )
        historyLen = kernel.tapsPer𝜙 - 1
    end

    history = zeros( historyLen )

    FIRFilter( kernel, history, historyLen )
end

function FIRFilter( h::Vector, rate::FloatingPoint, N𝜙::Integer = 32 )
    rate > 0.0 || error( "rate must be greater than 0" )
    kernel     = FIRArbitrary( h, rate, N𝜙 )
    historyLen = kernel.tapsPer𝜙 - 1
    history    = zeros( historyLen )
    FIRFilter( kernel, history, historyLen )
end




#==============================================================================#
#                            ____ ____ ____ ____ ___                           #
#                            |__/ |___ [__  |___  |                            #
#                            |  \ |___ ___] |___  |                            #
#==============================================================================#

# Resets filter and its kernel to an initial state

# Does nothing for non-rational kernels
reset( self::FIRKernel ) = self

# For rational kernel, set 𝜙Idx back to 1
reset( self::FIRRational ) = self.𝜙Idx = 1

# For rational kernel, set 𝜙Idx back to 1
function reset( self::FIRArbitrary )
    self.yCount = 0
    update!( self )
end

# For FIRFilter, set delay line to zeros of same tyoe and required length
function reset( self::FIRFilter )
    self.history = zeros( eltype( self.history ), self.historyLen )
    reset( self.kernel )
    return self
end




#==============================================================================#
#                      _  _    ___ ____    ___  ____ ___                       #
#                      |__|     |  |  |    |__] |___ |__]                      #
#                      |  |     |  |__|    |    |    |__]                      #
#==============================================================================#

# Converts a vector of coefficients to a matrix. Each column is a filter.
# Appends zeros if necessary.
# Example:
#   julia> taps2pfb( [1:9], 4 )
#   3x4 Array{Int64,2}:
#    1  2  3  4
#    5  6  7  8
#    9  0  0  0

function taps2pfb{T}( h::Vector{T}, N𝜙::Integer )
    hLen      = length( h )
    hLenPer𝜙  = iceil(  hLen/N𝜙  )
    pfbSize   = hLenPer𝜙 * N𝜙

    if hLen != pfbSize                                # check that the vector is an integer multiple of N𝜙
        hExtended             = similar( h, pfbSize ) # No? extend and zero pad
        hExtended[1:hLen]     = h
        hExtended[hLen+1:end] = 0
        h                     = hExtended
    end

    hLen      = length( h )
    hLenPer𝜙  = int( hLen/N𝜙 )
    pfb       = reshape( h, N𝜙, hLenPer𝜙 )'
end




#==============================================================================#
#               ____ _  _ ___ ___  _  _ ___    _    ____ _  _                  #
#               |  | |  |  |  |__] |  |  |     |    |___ |\ |                  #
#               |__| |__|  |  |    |__|  |     |___ |___ | \|                  #
#==============================================================================#

# Calculates the resulting length of a multirate filtering operation, given a
#   FIRFilter{FIRRational} object and an input vector
#
# ( It's hard to explain how this works without a diagram )

function outputlength( inputlength::Integer, ratio::Rational, initial𝜙::Integer )
    interpolation = num( ratio )
    decimation    = den( ratio )
    outLen        = (( inputlength * interpolation ) - initial𝜙 + 1 ) / decimation
    iceil(  outLen  )
end

function outputlength( self::FIRFilter{FIRStandard}, inputlength::Integer )
    inputlength
end

function outputlength( self::FIRFilter{FIRInterpolator}, inputlength::Integer )
    kernel = self.kernel
    kernel.interpolation * inputlength
end

function outputlength( self::FIRFilter{FIRDecimator}, inputlength::Integer )
    kernel = self.kernel
    outputlength( inputlength-kernel.inputDeficit+1, 1//kernel.decimation, 1 )
end

function outputlength( self::FIRFilter{FIRRational}, inputlength::Integer )
    kernel = self.kernel
    outputlength( inputlength-kernel.inputDeficit+1, kernel.ratio, kernel.𝜙Idx )
end

# TODO: outputlength function for arbitrary FIR kernel




#==============================================================================#
#                 _ _  _ ___  _  _ ___    _    ____ _  _                       #
#                 | |\ | |__] |  |  |     |    |___ |\ |                       #
#                 | | \| |    |__|  |     |___ |___ | \|                       #
#==============================================================================#

function inputlength( outputlength::Int, ratio::Rational, initial𝜙::Integer )
    interpolation = num( ratio )
    decimation    = den( ratio )
    inLen         = ( outputlength * decimation + initial𝜙 - 1 ) / interpolation
    iceil( inLen )
end

function inputlength( self::FIRFilter{FIRStandard}, outputlength::Integer )
    outputlength
end

function inputlength( self::FIRFilter{FIRInterpolator}, outputlength::Integer )
    kernel = self.kernel
    inputlength( outputlength, kernel.interpolation//1, 1 )
end

function inputlength( self::FIRFilter{FIRDecimator}, outputlength::Integer )
    kernel = self.kernel
    inLen  = inputlength( outputlength, 1//kernel.decimation, 1 )
    inLen  = inLen + kernel.inputlength - 1
end

function inputlength( self::FIRFilter{FIRRational}, outputlength::Integer )
    kernel = self.kernel
    inLen = inputlength( outputlength, kernel.ratio, kernel.𝜙Idx )
    inLen = inLen + kernel.inputDeficit - 1
end




#==============================================================================#
#              _  _ ____ _  _ ___    ___  _  _ ____ ____ ____                  #
#              |\ | |___  \/   |     |__] |__| |__| [__  |___                  #
#              | \| |___ _/\_  |     |    |  | |  | ___] |___                  #
#==============================================================================#

function nextphase( currentphase::Integer, ratio::Rational )
    interpolation = num( ratio )
    decimation    = den( ratio )
    𝜙Step         = mod( decimation, interpolation )
    𝜙Next         = currentphase + 𝜙Step
    𝜙Next         = 𝜙Next > interpolation ? 𝜙Next - interpolation : 𝜙Next
end




#==============================================================================#
#               ____ _ _  _ ____ _    ____    ____ ____ ___ ____               #
#               [__  | |\ | | __ |    |___    |__/ |__|  |  |___               #
#               ___] | | \| |__] |___ |___    |  \ |  |  |  |___               #
#==============================================================================#

function filt!{T}( buffer::Vector{T}, self::FIRFilter{FIRStandard}, x::Vector{T} )
    history::Vector{T} = self.history # TODO: figure out a better way of handling this
    h::Vector{T}       = self.kernel.h
    hLen               = self.kernel.hLen
    historyLen         = self.historyLen
    bufLen             = length( buffer )
    xLen               = length( x )
    outLen             = xLen
    criticalYidx       = min( hLen, outLen )

    bufLen >= xLen || error( "buffer length must be >= x length" )

    for yIdx in 1:criticalYidx        # this first loop takes care of filter ramp up and previous history
        @inbounds buffer[yIdx] = unsafedot( h, history, x, yIdx )
    end

    for yIdx in criticalYidx+1:xLen
        @inbounds buffer[yIdx] = unsafedot( h, x, yIdx )
    end

    self.history = shiftin!( history, x )

    return buffer
end

function filt{T}( self::FIRFilter{FIRStandard}, x::Vector{T} )
    buffer = zeros( eltype(x), length(x) )
    filt!( buffer, self, x )
end




#==============================================================================#
#               _ _  _ ___ ____ ____ ___  _    ____ ____ ___ ____              #
#               | |\ |  |  |___ |__/ |__] |    |  | |__|  |  |___              #
#               | | \|  |  |___ |  \ |    |___ |__| |  |  |  |___              #
#==============================================================================#

function filt!{T}( buffer::Vector{T}, self::FIRFilter{FIRInterpolator}, x::Vector{T} )
    pfb::PFB{T}        = self.kernel.pfb
    history::Vector{T} = self.history
    interpolation      = self.kernel.interpolation
    N𝜙                 = self.kernel.N𝜙
    tapsPer𝜙           = self.kernel.tapsPer𝜙
    xLen               = length( x )
    bufLen             = length( buffer )
    historyLen         = self.historyLen
    outLen             = outputlength( self, xLen )
    criticalYidx       = min( historyLen*interpolation, outLen )
    inputIdx           = 1
    𝜙                  = 1

    bufLen >= outLen || error( "length( buffer ) must be >= interpolation * length(x)")

    for yIdx in 1:criticalYidx
        @inbounds buffer[yIdx] = unsafedot( pfb, 𝜙, history, x, inputIdx )
        (𝜙, inputIdx) = 𝜙 == N𝜙 ? ( 1, inputIdx+1 ) : ( 𝜙+1, inputIdx )
    end
    for yIdx in criticalYidx+1:outLen
        @inbounds buffer[yIdx] = unsafedot( pfb, 𝜙, x, inputIdx )
        (𝜙, inputIdx) = 𝜙 == N𝜙 ? ( 1, inputIdx+1 ) : ( 𝜙+1, inputIdx )
    end

    self.history = shiftin!( history, x )

    return buffer
end

function filt( self::FIRFilter{FIRInterpolator}, x::Vector )
    xLen   = length( x )
    outlen = outputlength( self, xLen )
    buffer = similar( x, outlen )
    filt!( buffer, self, x )
end




#==============================================================================#
#           ____ ____ ___     ____ ____ ____ ____ _  _ ___  _    ____          #
#           |__/ |__|  |      |__/ |___ [__  |__| |\/| |__] |    |___          #
#           |  \ |  |  |  .   |  \ |___ ___] |  | |  | |    |___ |___          #
#==============================================================================#

function filt!{T}( buffer::Vector{T}, self::FIRFilter{FIRRational}, x::Vector{T} )
    kernel             = self.kernel
    xLen               = length( x )
    bufLen             = length( buffer )

    if xLen < kernel.inputDeficit
        self.history = shiftin!( history, x )
        kernel.inputDeficit -= xLen
        return T[]
    end

    outLen = outputlength( xLen-kernel.inputDeficit+1, kernel.ratio, kernel.𝜙Idx )
    bufLen >= outLen || error( "buffer is too small" )

    pfb::PFB{T}        = kernel.pfb
    history::Vector{T} = self.history
    interpolation      = num( kernel.ratio )
    decimation         = den( kernel.ratio )
    𝜙IdxStepSize       = mod( decimation, interpolation )
    critical𝜙Idx       = kernel.N𝜙 - 𝜙IdxStepSize
    inputIdx           = kernel.inputDeficit
    yIdx               = 0

    while inputIdx <= xLen
        yIdx += 1
        if inputIdx < kernel.tapsPer𝜙
            hIdx = 1
            accumulator = unsafedot( pfb, kernel.𝜙Idx, history, x, inputIdx )
        else
            accumulator = unsafedot( pfb, kernel.𝜙Idx, x, inputIdx )
        end

        buffer[ yIdx ] = accumulator
        inputIdx      += ifloor( ( kernel.𝜙Idx + decimation - 1 ) / interpolation )
        kernel.𝜙Idx    = nextphase( kernel.𝜙Idx, kernel.ratio )
    end

    kernel.inputDeficit = inputIdx - xLen
    self.history        = shiftin!( history, x )

    return yIdx
end

function filt{T}( self::FIRFilter{FIRRational}, x::Vector{T} )
    kernel = self.kernel
    xLen   = length( x )

    if xLen < kernel.inputDeficit
        history::Vector{T} = self.history
        self.history = shiftin!( history, x )
        kernel.inputDeficit -= xLen
        return T[]
    end

    outLen = outputlength( self, xLen )
    buffer = similar( x, outLen )
    filt!( buffer, self, x )

    return buffer
end




#==============================================================================#
#                      ___  ____ ____ _ _  _ ____ ___ ____                     #
#                      |  \ |___ |    | |\/| |__|  |  |___                     #
#                      |__/ |___ |___ | |  | |  |  |  |___                     #
#==============================================================================#

function filt!{T}( buffer::Vector{T}, self::FIRFilter{FIRDecimator}, x::Vector{T} )
    kernel = self.kernel
    xLen   = length( x )

    if xLen < kernel.inputDeficit
        self.history = shiftin!( history, x )
        kernel.inputDeficit -= xLen
        return T[]
    end

    outLen = outputlength( self, xLen )
    history::Vector{T} = self.history
    h::Vector{T}       = kernel.h
    inputIdx           = kernel.inputDeficit
    yIdx               = 0

    while inputIdx <= xLen
        accumulator = zero( T )
        yIdx       += 1

        if inputIdx < kernel.hLen
            accumulator = unsafedot( h, history, x, inputIdx )
        else
            accumulator = unsafedot( h, x, inputIdx )
        end

        buffer[ yIdx ] = accumulator
        inputIdx      += kernel.decimation
    end

    kernel.inputDeficit = inputIdx - xLen
    self.history        = shiftin!( history, x )

    return yIdx
end

function filt{T}( self::FIRFilter{FIRDecimator}, x::Vector{T} )
    kernel             = self.kernel
    xLen               = length( x )

    if xLen < kernel.inputDeficit
        history::Vector{T} = self.history
        self.history       = shiftin!( history, x )
        kernel.inputDeficit -= xLen
        return T[]
    end

    outLen = outputlength( self, xLen )
    buffer = similar( x, outLen )
    filt!( buffer, self, x )

    return buffer
end




#==============================================================================#
#        ____ ____ ___      ____ ____ ____ ____ _  _ ___  _    ____ ____       #
#        |__| |__/ |__]     |__/ |___ [__  |__| |\/| |__] |    |___ |__/       #
#        |  | |  \ |__] .   |  \ |___ ___] |  | |  | |    |___ |___ |  \       #
#==============================================================================#

# Updates FIRArbitrary state. See Section 7.5.1 in [1].
#   [1] uses a phase accumilator that increments by Δ (N𝜙/rate)
#   The original implementation of update! used this method, but the numerical
#   errors built up pretty quickly. Instead α is now incremented by δ, the
#   fractional part of Δ. The phase index, 𝜙Idx, is now incremented by
#   integer part of Δ plus the integer part of α. However, α should always
#   be < 1. So α is first incremented by δ and may be greater than 1 at this
#   point. We dont want to fix that until we add the integer part to 𝜙Idx first.
#   I hope that makes sense.
function update!( kernel::FIRArbitrary )
    kernel.α    += kernel.δ
    kernel.𝜙Idx += ifloor( kernel.α ) + kernel.𝜙Stride
    kernel.α     = mod( kernel.α, 1.0 )

    if kernel.𝜙Idx > kernel.N𝜙
        kernel.xIdx += ifloor( (kernel.𝜙Idx-1) / kernel.N𝜙 )
        kernel.𝜙Idx  = mod( (kernel.𝜙Idx-1), kernel.N𝜙 ) + 1
    end
end

function filt{Th,Tx}( self::FIRFilter{FIRArbitrary{Th}}, x::Vector{Tx} )
    kernel              = self.kernel
    pfb                 = kernel.pfb
    dpfb                = kernel.dpfb
    xLen                = length( x )
    bufLen              = iceil( xLen * kernel.rate ) + 1
    buffer              = Array(promote_type(Th,Tx), bufLen)
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
    # We do this by seting inputIdx to inputDeficit which was calculated in the previous run.
    # InputDeficit is set to 1 when instantiation the FIRArbitrary kernel, that way the first
    #   input always produces an output.
    kernel.xIdx = kernel.inputDeficit

    while kernel.xIdx <= xLen
        db_vec_xidx[bufIdx] = kernel.xIdx
        db_vec_phi[bufIdx]  = kernel.𝜙Idx + kernel.α
        if kernel.xIdx < kernel.tapsPer𝜙
            yLower = unsafedot( pfb,  kernel.𝜙Idx, history, x, kernel.xIdx )
            yUpper = unsafedot( dpfb, kernel.𝜙Idx, history, x, kernel.xIdx )
        else
            yLower = unsafedot( pfb,  kernel.𝜙Idx, x, kernel.xIdx )
            yUpper = unsafedot( dpfb, kernel.𝜙Idx, x, kernel.xIdx )
        end
        buffer[bufIdx] = yLower + yUpper * kernel.α
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
    # return buffer
end




#==============================================================================#
#              ____ ____ ____ ____ ____ _ _ _    ____ _ _    ___               #
#              |___ |__| |__/ |__/ |  | | | |    |___ | |     |                #
#              |    |  | |  \ |  \ |__| |_|_|    |    | |___  |                #
#==============================================================================#

include( "FIRFarrow.jl" )




#==============================================================================#
#       ____ ___ ____ ___ ____ _    ____ ____ ____    ____ _ _    ___          #
#       [__   |  |__|  |  |___ |    |___ [__  [__     |___ | |     |           #
#       ___]  |  |  |  |  |___ |___ |___ ___] ___]    |    | |___  |           #
#==============================================================================#

function filt( h::Vector, x::Vector, ratio::Rational = 1//1 )
    self = FIRFilter( h, ratio )
    filt( self, x )
end

function filt( h::Vector, x::Vector, rate::FloatingPoint, N𝜙::Integer = 32 )
    self = FIRFilter( h, rate, N𝜙 )
    filt( self, x )
end




#==============================================================================#
#                   ____ _ ___ ____ ___ _ ____ _  _ ____                       #
#                   |    |  |  |__|  |  | |  | |\ | [__                        #
#                   |___ |  |  |  |  |  | |__| | \| ___]                       #
#==============================================================================#

# [1] F.J. Harris, *Multirate Signal Processing for Communication Systems*. Prentice Hall, 2004
# [2] Dick, C.; Harris, F., "Options for arbitrary resamplers in FPGA-based modulators," Signals, Systems and Computers, 2004. Conference Record of the Thirty-Eighth Asilomar Conference on , vol.1, no., pp.777,781 Vol.1, 7-10 Nov. 2004
# [3] Kim, S.C.; Plishker, W.L.; Bhattacharyya, S.S., "An efficient GPU implementation of an arbitrary resampling polyphase channelizer," Design and Architectures for Signal and Image Processing (DASIP), 2013 Conference on, vol., no., pp.231,238, 8-10 Oct. 2013
