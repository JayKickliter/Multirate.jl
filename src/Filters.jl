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
end

function FIRStandard( h::Vector )
    h    = flipud( h )
    hLen = length( h )
    FIRStandard( h, hLen )
end


# Interpolator FIR kernel
type FIRInterpolator <: FIRKernel
    pfb::PFB
    interpolation::Int
    N𝜙::Int
    tapsPer𝜙::Int
end

function FIRInterpolator( h::Vector, interpolation::Integer )
    pfb              = flipud( polyize( h, interpolation ) )
    ( tapsPer𝜙, N𝜙 ) = size( pfb )
    FIRInterpolator( pfb, interpolation, N𝜙, tapsPer𝜙 )
end


# Decimator FIR kernel
type FIRDecimator <: FIRKernel
    h::Vector
    hLen::Int
    decimation::Int
    inputDeficit::Int
end

function FIRDecimator( h::Vector, decimation::Integer )
    h    = flipud( h )
    hLen = length( h )
    FIRDecimator( h, hLen, decimation, 1 )
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
end

function FIRRational( h::Vector, ratio::Rational )
    interpolation    = num( ratio )
    pfb              = flipud( polyize( h, interpolation ) )
    ( tapsPer𝜙, N𝜙 ) = size( pfb )
    criticalYidx     = ifloor( tapsPer𝜙 * ratio )
    FIRRational( pfb, ratio, N𝜙, tapsPer𝜙, criticalYidx, 1, 1 )
end


# Arbitrary resampler FIR kernel
type FIRArbitrary  <: FIRKernel
    pfb::PFB
    N𝜙::Int
    tapsPer𝜙::Int
    resampleRate::Float64
    yCount::Int
    xCount::Int
    yLower::Number
    yUpperStalled::Bool
    𝜙IdxLower::Int
    𝜙IdxUpper::Int
    xIdxDelta::Int
    xIdxUpperOffset::Int
    inputDeficit::Int
    Δ::Float64
    ΔPrevious::Float64
    function FIRArbitrary( h::Vector, resampleRate::Real, numFilters::Integer )
        pfb             = flipud( polyize( h, numFilters ) )
        tapsPer𝜙        = size( pfb )[1]
        N𝜙              = size( pfb )[2]
        resampleRate    = resampleRate
        yCount          = 0
        xCount          = 0
        yLower          = NaN
        yUpperStalled   = false
        𝜙IdxLower       = 0
        𝜙IdxUpper       = 0
        inputDeficit    = 1
        xIdxDelta       = 0
        xIdxUpperOffset = 0
        Δ               = 0.0
        ΔPrevious       = 0.0
        new( pfb, N𝜙, tapsPer𝜙, resampleRate, yCount, xCount, yLower, yUpperStalled, 𝜙IdxLower, 𝜙IdxUpper, xIdxDelta, xIdxUpperOffset, inputDeficit, Δ, ΔPrevious )
    end
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
    historyLen = 0

    if resampleRatio == 1                                     # single-rate
        kernel        = FIRStandard( h )
        historyLen = kernel.hLen - 1
    elseif interpolation == 1                                 # decimate
        kernel        = FIRDecimator( h, decimation )
        historyLen = kernel.hLen - 1
    elseif decimation == 1                                    # interpolate
        kernel        = FIRInterpolator( h, interpolation )
        historyLen = kernel.tapsPer𝜙 - 1
    else                                                      # rational
        kernel        = FIRRational( h, resampleRatio )
        historyLen = kernel.tapsPer𝜙 - 1
    end

    history = zeros( historyLen )

    FIRFilter( kernel, history, historyLen )
end

function FIRFilter( h::Vector, resampleRate::FloatingPoint, numFilters::Integer = 32 )
    resampleRate > 0.0 || error( "resampleRate must be greater than 0" )

    kernel        = FIRArbitrary( h, resampleRate, numFilters )
    historyLen = kernel.tapsPer𝜙 - 1
    history       = zeros( historyLen )

    updatestate!( kernel )

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
#   julia> polyize( [1:9], 4 )
#   3x4 Array{Int64,2}:
#    1  2  3  4
#    5  6  7  8
#    9  0  0  0

function polyize{T}( h::Vector{T}, numFilters::Integer )
    hLen      = length( h )
    hLenPer𝜙  = iceil(  hLen/numFilters  )
    pfbSize   = hLenPer𝜙 * numFilters

    if hLen != pfbSize                                # check that the vector is an integer multiple of numFilters
        hExtended             = similar( h, pfbSize ) # No? extend and zero pad
        hExtended[1:hLen]     = h
        hExtended[hLen+1:end] = 0
        h                     = hExtended
    end

    hLen      = length( h )
    hLenPer𝜙  = int( hLen/numFilters )
    pfb       = reshape( h, numFilters, hLenPer𝜙 )'
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




#==============================================================================#
#                 _ _  _ ___  _  _ ___    _    ____ _  _                       #
#                 | |\ | |__] |  |  |     |    |___ |\ |                       #
#                 | | \| |    |__|  |     |___ |___ | \|                       #
#==============================================================================#

function inputlength( outputlength::Int, ratio::Rational, initial𝜙::Integer )
    interpolation = num( ratio )
    decimation    = den( ratio )
    inLen = ( outputlength * decimation + initial𝜙 - 1 ) / interpolation
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
    history::Vector{T} = self.history
    h::Vector{T}       = self.kernel.h
    hLen               = self.kernel.hLen
    historyLen         = self.historyLen
    bufLen             = length( buffer )
    xLen               = length( x )
    outLen             = xLen
    criticalYidx       = min( hLen, outLen )

    bufLen >= xLen || error( "buffer length must be >= x length" )

    for yIdx in 1:criticalYidx                                   # this first loop takes care of filter ramp up and previous history
        @inbounds buffer[yIdx] = unsafedot( h, history, x, yIdx )
    end

    for yIdx in criticalYidx+1:xLen
        @inbounds buffer[yIdx] = unsafedot( h, x, yIdx )
    end

    if xLen >= self.historyLen
        copy!( history, 1, x, xLen - self.historyLen + 1, self.historyLen )
    else
        history = [ history, x ][ end - self.historyLen + 1: end ]
    end
    self.history = history

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

    bufLen >= outLen || error( "length( buffer ) must be >= interpolation * length(x)")

    inputIdx = 1
    𝜙        = 1

    for yIdx in 1:criticalYidx
        @inbounds buffer[yIdx] = unsafedot( pfb, 𝜙, history, x, inputIdx )
        (𝜙, inputIdx) = 𝜙 == N𝜙 ? ( 1, inputIdx+1 ) : ( 𝜙+1, inputIdx )
    end

    for yIdx in criticalYidx+1:outLen
        @inbounds buffer[yIdx] = unsafedot( pfb, 𝜙, x, inputIdx )
        (𝜙, inputIdx) = 𝜙 == N𝜙 ? ( 1, inputIdx+1 ) : ( 𝜙+1, inputIdx )
    end

    if xLen >= self.historyLen
        copy!( history, 1, x, xLen - self.historyLen + 1, self.historyLen )
    else
        history = [ history, x ][ end - self.historyLen + 1: end ]
    end
    self.history = history


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
    kernel = self.kernel
    xLen   = length( x )
    bufLen = length( buffer )

    if xLen < kernel.inputDeficit
        self.history = [ self.history, x ][ end - self.historyLen + 1: end ]
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

        inputIdx   += ifloor( ( kernel.𝜙Idx + decimation - 1 ) / interpolation )
        kernel.𝜙Idx = nextphase( kernel.𝜙Idx, kernel.ratio )
    end

    kernel.inputDeficit = inputIdx - xLen

    if xLen >= self.historyLen
        copy!( history, 1, x, xLen - self.historyLen + 1, self.historyLen )
    else
        history = [ history, x ][ end - self.historyLen + 1: end ]
    end
    self.history = history

    return yIdx
end

function filt{T}( self::FIRFilter{FIRRational}, x::Vector{T} )
    kernel = self.kernel
    xLen   = length( x )

    if xLen < kernel.inputDeficit
        self.history = [ self.history, x ][ end - self.historyLen + 1: end ]
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
        self.history = [ self.history, x ][ end - self.historyLen + 1: end ]
        kernel.inputDeficit -= xLen
        return T[]
    end

    outLen = outputlength( self, xLen )

    h::Vector{T}       = kernel.h
    history::Vector{T} = self.history
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

    if xLen >= self.historyLen
        copy!( history, 1, x, xLen - self.historyLen + 1, self.historyLen )
    else
        history = [ history, x ][ end - self.historyLen + 1: end ]
    end
    self.history = history


    return yIdx
end

function filt{T}( self::FIRFilter{FIRDecimator}, x::Vector{T} )
    kernel = self.kernel
    xLen   = length( x )

    if xLen < kernel.inputDeficit
        self.history = [ self.history, x ][ end - self.historyLen + 1: end ]
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

function updatestate!( self::FIRArbitrary )
    self.yCount         += 1
    self.𝜙IdxLower       = ifloor( mod( (self.yCount-1)/self.resampleRate, 1 ) * self.N𝜙 ) + 1
    self.𝜙IdxUpper       = self.𝜙IdxLower == self.N𝜙 ? 1 : self.𝜙IdxLower + 1
    self.xIdxUpperOffset = self.𝜙IdxLower == self.N𝜙 ? 1 : 0
    xCountCurrent        = self.xCount
    self.xCount          = ifloor( (self.yCount-1)/self.resampleRate )
    self.xIdxDelta       = self.xCount - xCountCurrent
    self.ΔPrevious       = self.Δ
    self.Δ               = mod( (self.yCount-1) * self.N𝜙 / self.resampleRate, 1 )
end

function filt{T}( self::FIRFilter{FIRArbitrary}, x::Vector{T} )
    kernel             = self.kernel
    xLen               = length( x )
    bufLen             = iceil( xLen * kernel.resampleRate ) + 1
    buffer             = similar( x, bufLen )
    bufIdx             = 1
    history::Vector{T} = self.history
    pfb::PFB{T}        = kernel.pfb

    # In the previous run, did 𝜙IdxUpper wrap around, requiring an extra input that we didn't have yet?
    if kernel.yUpperStalled && xLen >= 1
        yUpper               = dot( pfb[:,1], [ self.history, x[1] ]  )
        buffer[bufIdx]       = kernel.yLower * (1 - kernel.ΔPrevious) + yUpper * kernel.ΔPrevious
        kernel.yUpperStalled = false
        bufIdx              += 1
    end

    # Do we have enough input samples to produce one or more output samples?
    if xLen < kernel.inputDeficit
        self.history = [ self.history, x ][ end - self.historyLen + 1: end ]
        kernel.inputDeficit -= xLen
        return buffer[1:bufIdx-1]
    end

    # Skip over input samples that are not needed to produce output results.
    # We do this by seting inputIdx to inputDeficit which was calculated in the previous run.
    # InputDeficit is set to 1 when instantiation the FIRArbitrary kernel, that way the first
    #   input always produces an output.
    inputIdx = kernel.inputDeficit

    while inputIdx <= xLen
        yLower = zero(T)
        yUpper = zero(T)
        xIdx   = inputIdx

        # Compute yLower
        #   As long as inputIdx <= xLen, we can calculate yLower
        if xIdx < kernel.tapsPer𝜙
            yLower = unsafedot( pfb, kernel.𝜙IdxLower, history, x, xIdx )
        else
            yLower = unsafedot( pfb, kernel.𝜙IdxLower, x, xIdx )
        end

        # If 𝜙IdxUpper wraps around, we will need another input sample.
        kernel.yLower = yLower
        xIdx += kernel.xIdxUpperOffset

        # If xIdx was advanced in the previous line due to a 𝜙IdxUpper wrap-around,
        #   we need to make sure we there are still enough input samples to complete this output.
        if xIdx <= xLen
            if xIdx < kernel.tapsPer𝜙
                yUpper = unsafedot( pfb, kernel.𝜙IdxUpper, history, x, xIdx )
            else
                yUpper = unsafedot( pfb, kernel.𝜙IdxUpper, x, xIdx )
            end
            buffer[bufIdx] = yLower * (1 - kernel.Δ) + yUpper * kernel.Δ
            bufIdx   += 1
        else
            # To finish computing this output sample, we need to compute yUpper.
            # However, we've reached the end of the line.
            # Set the 'stalled' state in the kernel and finish this output next time.
            kernel.yUpperStalled = true
        end
        updatestate!( kernel )
        inputIdx += kernel.xIdxDelta
    end

    # Did we overestimate needed buffer size?
    # TODO: Get rid of this by correctly calculating output size.
    bufLen == bufIdx - 1 || resize!( buffer, bufIdx - 1)
    kernel.inputDeficit = inputIdx - xLen

    if xLen >= self.historyLen
        copy!( history, 1, x, xLen - self.historyLen + 1, self.historyLen )
    else
        history = [ history, x ][ end - self.historyLen + 1: end ]
    end
    self.history = history

    return buffer
end




#==============================================================================#
#       ____ ___ ____ ___ ____ _    ____ ____ ____    ____ _ _    ___          #
#       [__   |  |__|  |  |___ |    |___ [__  [__     |___ | |     |           #
#       ___]  |  |  |  |  |___ |___ |___ ___] ___]    |    | |___  |           #
#==============================================================================#

function filt( h::Vector, x::Vector, ratio::Rational = 1//1 )
    self = FIRFilter( h, ratio )
    filt( self, x )
end
