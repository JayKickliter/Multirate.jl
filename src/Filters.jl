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
    α::Float64
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
        α               = 0.0
        new( pfb, N𝜙, tapsPer𝜙, resampleRate, yCount, xCount, yLower, yUpperStalled, 𝜙IdxLower, 𝜙IdxUpper, xIdxDelta, xIdxUpperOffset, inputDeficit, α )
    end
end



# FIRFilter - the kernel does the heavy lifting
type FIRFilter{Tk<:FIRKernel} <: Filter
    kernel::Tk
    dlyLine::Vector
    reqDlyLineLen::Int
end

function FIRFilter( h::Vector, resampleRatio::Rational = 1//1 )
    interpolation = num( resampleRatio )
    decimation    = den( resampleRatio )
    reqDlyLineLen = 0

    if resampleRatio == 1                                     # single-rate
        kernel        = FIRStandard( h )
        reqDlyLineLen = kernel.hLen - 1
    elseif interpolation == 1                                 # decimate
        kernel        = FIRDecimator( h, decimation )
        reqDlyLineLen = kernel.hLen - 1
    elseif decimation == 1                                    # interpolate
        kernel        = FIRInterpolator( h, interpolation )
        reqDlyLineLen = kernel.tapsPer𝜙 - 1
    else                                                      # rational
        kernel        = FIRRational( h, resampleRatio )
        reqDlyLineLen = kernel.tapsPer𝜙 - 1
    end

    dlyLine = zeros( reqDlyLineLen )

    FIRFilter( kernel, dlyLine, reqDlyLineLen )
end

function FIRFilter( h::Vector, resampleRate::FloatingPoint, numFilters::Integer = 32 )
    resampleRate > 0.0 || error( "resampleRate must be greater than 0" )

    kernel        = FIRArbitrary( h, resampleRate, numFilters )
    reqDlyLineLen = kernel.tapsPer𝜙 - 1
    dlyLine       = zeros( reqDlyLineLen )

    updatestate!( kernel )

    FIRFilter( kernel, dlyLine, reqDlyLineLen )
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
    self.dlyLine = zeros( eltype( self.dlyLine ), self.reqDlyLineLen )
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
    dlyLine::Vector{T} = self.dlyLine
    h::Vector{T}       = self.kernel.h
    hLen               = self.kernel.hLen
    reqDlyLineLen      = self.reqDlyLineLen
    bufLen             = length( buffer )
    xLen               = length( x )
    outLen             = xLen
    criticalYidx       = min( hLen, outLen )

    bufLen >= xLen || error( "buffer length must be >= x length" )

    for yIdx in 1:criticalYidx                                   # this first loop takes care of filter ramp up and previous dlyLine
        accumulator = zero(T)

        for k in 1:hLen-yIdx
            @inbounds accumulator += h[k] * dlyLine[k+yIdx-1]
        end

        for k in 1:yIdx
            @inbounds accumulator += h[hLen-yIdx+k] * x[k]
        end

        @inbounds buffer[yIdx] = accumulator
    end

    for yIdx in criticalYidx+1:xLen
        accumulator = zero(T)

        for k in 1:hLen
            @inbounds accumulator += h[k] * x[yIdx-hLen+k]
        end

        @inbounds buffer[yIdx] = accumulator
    end

    if xLen >= self.reqDlyLineLen
        copy!( dlyLine, 1, x, xLen - self.reqDlyLineLen + 1, self.reqDlyLineLen )
    else
        dlyLine = [ dlyLine, x ][ end - self.reqDlyLineLen + 1: end ]
    end
    self.dlyLine = dlyLine

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
    dlyLine::Vector{T} = self.dlyLine
    interpolation      = self.kernel.interpolation
    N𝜙                 = self.kernel.N𝜙
    tapsPer𝜙           = self.kernel.tapsPer𝜙
    xLen               = length( x )
    bufLen             = length( buffer )
    reqDlyLineLen      = self.reqDlyLineLen
    outLen             = outputlength( self, xLen )
    criticalYidx       = min( reqDlyLineLen*interpolation, outLen )

    bufLen >= outLen || error( "length( buffer ) must be >= interpolation * length(x)")

    inputIdx = 1
    𝜙        = 1

    for yIdx in 1:criticalYidx

        accumulator = zero(T)

        for k in 1:tapsPer𝜙-inputIdx
            @inbounds accumulator += pfb[k, 𝜙] * dlyLine[k+inputIdx-1]
        end

        for k in 1:inputIdx
            @inbounds accumulator += pfb[tapsPer𝜙-inputIdx+k, 𝜙] * x[k]
        end

        @inbounds buffer[yIdx]  = accumulator
        (𝜙, inputIdx) = 𝜙 == N𝜙 ? ( 1, inputIdx+1 ) : ( 𝜙+1, inputIdx )
    end

    for yIdx in criticalYidx+1:outLen

        accumulator = zero(T)

        for k in 1:tapsPer𝜙
            @inbounds accumulator += pfb[ k, 𝜙 ] * x[ inputIdx - tapsPer𝜙 + k ]
        end

        @inbounds buffer[yIdx]  = accumulator
        (𝜙, inputIdx) = 𝜙 == N𝜙 ? ( 1, inputIdx+1 ) : ( 𝜙+1, inputIdx )
    end

    if xLen >= self.reqDlyLineLen
        copy!( dlyLine, 1, x, xLen - self.reqDlyLineLen + 1, self.reqDlyLineLen )
    else
        dlyLine = [ dlyLine, x ][ end - self.reqDlyLineLen + 1: end ]
    end
    self.dlyLine = dlyLine


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
        self.dlyLine = [ self.dlyLine, x ][ end - self.reqDlyLineLen + 1: end ]
        kernel.inputDeficit -= xLen
        return T[]
    end

    outLen = outputlength( xLen-kernel.inputDeficit+1, kernel.ratio, kernel.𝜙Idx )
    bufLen >= outLen || error( "buffer is too small" )

    pfb::PFB{T}        = kernel.pfb
    dlyLine::Vector{T} = self.dlyLine
    interpolation      = num( kernel.ratio )
    decimation         = den( kernel.ratio )
    𝜙IdxStepSize       = mod( decimation, interpolation )
    critical𝜙Idx       = kernel.N𝜙 - 𝜙IdxStepSize

    inputIdx           = kernel.inputDeficit
    yIdx               = 0

    while inputIdx <= xLen

        accumulator = zero( T )
        yIdx       += 1

        if inputIdx < kernel.tapsPer𝜙
            hIdx = 1
            for k in inputIdx:self.reqDlyLineLen
                @inbounds accumulator += pfb[ hIdx, kernel.𝜙Idx ] * dlyLine[ k ]
                hIdx += 1
            end

            for k in 1:inputIdx
                @inbounds accumulator += pfb[ hIdx, kernel.𝜙Idx ] * x[ k ]
                hIdx += 1
            end
        else
            hIdx = 1
            for k in inputIdx-kernel.tapsPer𝜙+1:inputIdx
                @inbounds accumulator += pfb[ hIdx, kernel.𝜙Idx ] * x[ k ]
                hIdx += 1
            end
        end

        buffer[ yIdx ] = accumulator

        inputIdx   += ifloor( ( kernel.𝜙Idx + decimation - 1 ) / interpolation )
        kernel.𝜙Idx = nextphase( kernel.𝜙Idx, kernel.ratio )
    end

    kernel.inputDeficit = inputIdx - xLen

    if xLen >= self.reqDlyLineLen
        copy!( dlyLine, 1, x, xLen - self.reqDlyLineLen + 1, self.reqDlyLineLen )
    else
        dlyLine = [ dlyLine, x ][ end - self.reqDlyLineLen + 1: end ]
    end
    self.dlyLine = dlyLine

    return yIdx
end

function filt{T}( self::FIRFilter{FIRRational}, x::Vector{T} )
    kernel = self.kernel
    xLen   = length( x )

    if xLen < kernel.inputDeficit
        self.dlyLine = [ self.dlyLine, x ][ end - self.reqDlyLineLen + 1: end ]
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
        self.dlyLine = [ self.dlyLine, x ][ end - self.reqDlyLineLen + 1: end ]
        kernel.inputDeficit -= xLen
        return T[]
    end

    outLen = outputlength( self, xLen )

    h::Vector{T}       = kernel.h
    dlyLine::Vector{T} = self.dlyLine
    inputIdx           = kernel.inputDeficit
    yIdx               = 0

    while inputIdx <= xLen

        accumulator = zero( T )
        yIdx       += 1

        if inputIdx < kernel.hLen
            hIdx = 1
            for k in inputIdx:self.reqDlyLineLen
                @inbounds accumulator += h[ hIdx ] * dlyLine[ k ]
                hIdx += 1
            end

            for k in 1:inputIdx
                @inbounds accumulator += h[ hIdx ] * x[ k ]
                hIdx += 1
            end
        else
            hIdx = 1
            for k in inputIdx-kernel.hLen+1:inputIdx
                @inbounds accumulator += h[ hIdx ] * x[ k ]
                hIdx += 1
            end
        end

        buffer[ yIdx ] = accumulator

        inputIdx   += kernel.decimation
    end

    kernel.inputDeficit = inputIdx - xLen

    if xLen >= self.reqDlyLineLen
        copy!( dlyLine, 1, x, xLen - self.reqDlyLineLen + 1, self.reqDlyLineLen )
    else
        dlyLine = [ dlyLine, x ][ end - self.reqDlyLineLen + 1: end ]
    end
    self.dlyLine = dlyLine


    return yIdx
end

function filt{T}( self::FIRFilter{FIRDecimator}, x::Vector{T} )
    kernel = self.kernel
    xLen   = length( x )

    if xLen < kernel.inputDeficit
        self.dlyLine = [ self.dlyLine, x ][ end - self.reqDlyLineLen + 1: end ]
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
    self.α               = mod( (self.yCount-1) * self.N𝜙 / self.resampleRate, 1 )
end

function filt{T}( self::FIRFilter{FIRArbitrary}, x::Vector{T} )
    kernel             = self.kernel
    xLen               = length( x )
    bufLen             = iceil( xLen * kernel.resampleRate ) + 1
    buffer             = similar( x, bufLen )
    pfb::PFB{T}        = kernel.pfb
    dlyLine::Vector{T} = self.dlyLine
    bufIdx               = 1

    if kernel.yUpperStalled && xLen >= 1
        yUpper       = dot( kernel.pfb[:,1], [ self.dlyLine, x[1] ]  )
        buffer[bufIdx] = kernel.yLower * (1 - kernel.α) + yUpper * kernel.α
        bufIdx        += 1
    end

    kernel.yUpperStalled

    if xLen < kernel.inputDeficit
        self.dlyLine = [ self.dlyLine, x ][ end - self.reqDlyLineLen + 1: end ]
        kernel.inputDeficit -= xLen
        return buffer
    end

    inputIdx = kernel.inputDeficit

    while inputIdx <= xLen
        println( "yCount = $(kernel.yCount), 𝜙Idx = $(kernel.𝜙IdxLower), inputIdx = $inputIdx, α = $(kernel.α)" )
        yLower = zero(T)
        yUpper = zero(T)

        xIdx = inputIdx

        # Compute yLower
        if xIdx < kernel.tapsPer𝜙
            hIdx = 1
            for k in xIdx:self.reqDlyLineLen
                yLower += pfb[ hIdx, kernel.𝜙IdxLower ] * dlyLine[ k ]
                hIdx += 1
            end
            for k in 1:xIdx
                yLower += pfb[ hIdx, kernel.𝜙IdxLower ] * x[ k ]
                hIdx += 1
            end
        else
            hIdx = 1
            for k in xIdx-kernel.tapsPer𝜙+1:xIdx
                yLower += pfb[ hIdx, kernel.𝜙IdxLower ] * x[ k ]
                hIdx += 1
            end
        end

        xIdx += kernel.xIdxUpperOffset
        if xIdx <= xLen
            # Compute yUpper

            if xIdx < kernel.tapsPer𝜙
                hIdx = 1
                for k in xIdx:self.reqDlyLineLen
                    yUpper += pfb[ hIdx, kernel.𝜙IdxUpper ] * dlyLine[ k ]
                    hIdx += 1
                end
                for k in 1:xIdx
                    yUpper += pfb[ hIdx, kernel.𝜙IdxUpper ] * x[ k ]
                    hIdx += 1
                end
            else
                hIdx = 1
                for k in xIdx-kernel.tapsPer𝜙+1:xIdx
                    yUpper += pfb[ hIdx, kernel.𝜙IdxUpper ] * x[ k ]
                    hIdx += 1
                end
            end

            buffer[bufIdx] = yLower * (1 - kernel.α) + yUpper * kernel.α
        else
            kernel.yLower        = yLower
            kernel.yUpperStalled = true
        end

        updatestate!( kernel )
        inputIdx += kernel.xIdxDelta
        bufIdx   += 1
    end

    resize!( buffer, bufIdx - 1)

    kernel.inputDeficit = inputIdx - xLen

    if xLen >= self.reqDlyLineLen
        copy!( dlyLine, 1, x, xLen - self.reqDlyLineLen + 1, self.reqDlyLineLen )
    else
        dlyLine = [ dlyLine, x ][ end - self.reqDlyLineLen + 1: end ]
    end
    self.dlyLine = dlyLine

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
