#==============================================================================#
#       ___ _   _ ___  ____ ____      /    ____ ____ _  _ ____ ___ ____        #
#        |   \_/  |__] |___ [__      /     |    |  | |\ | [__   |  |__/        #
#        |    |   |    |___ ___]    /      |___ |__| | \| ___]  |  |  \ .      #
#                                                                              #
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
    Nφ::Int
    tapsPerφ::Int
end

function FIRInterpolator( h::Vector, interpolation::Integer )
    pfb              = flipud( polyize( h, interpolation ) )
    ( tapsPerφ, Nφ ) = size( pfb )
    FIRInterpolator( pfb, interpolation, Nφ, tapsPerφ )
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
    Nφ::Int
    tapsPerφ::Int
    criticalYidx::Int
    φIdx::Int
    inputDeficit::Int
end

function FIRRational( h::Vector, ratio::Rational )
    interpolation    = num( ratio )
    pfb              = flipud( polyize( h, interpolation ) )
    ( tapsPerφ, Nφ ) = size( pfb )
    criticalYidx     = ifloor( tapsPerφ * ratio )
    FIRRational( pfb, ratio, Nφ, tapsPerφ, criticalYidx, 1, 1 )
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
        reqDlyLineLen = kernel.tapsPerφ - 1
    else                                                      # rational
        kernel        = FIRRational( h, resampleRatio )
        reqDlyLineLen = kernel.tapsPerφ - 1
    end

    dlyLine = zeros( reqDlyLineLen )

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

# For rational kernel, set φIdx back to 1
reset( self::FIRRational ) = self.φIdx = 1

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
    hLenPerφ  = iceil(  hLen/numFilters  )
    pfbSize   = hLenPerφ * numFilters

    if hLen != pfbSize                                # check that the vector is an integer multiple of numFilters
        hExtended             = similar( h, pfbSize ) # No? extend and zero pad
        hExtended[1:hLen]     = h
        hExtended[hLen+1:end] = 0
        h                     = hExtended
    end

    hLen      = length( h )
    hLenPerφ  = int( hLen/numFilters )
    pfb       = reshape( h, numFilters, hLenPerφ )'
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

function outputlength( inputlength::Integer, ratio::Rational, initialφ::Integer )
    interpolation = num( ratio )
    decimation    = den( ratio )
    outLen        = (( inputlength * interpolation ) - initialφ + 1 ) / decimation
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
    outputlength( inputlength-kernel.inputDeficit+1, kernel.ratio, kernel.φIdx )
end




#==============================================================================#
#                 _ _  _ ___  _  _ ___    _    ____ _  _                       #
#                 | |\ | |__] |  |  |     |    |___ |\ |                       #
#                 | | \| |    |__|  |     |___ |___ | \|                       #
#==============================================================================#

function inputlength( outputlength::Int, ratio::Rational, initialφ::Integer )
    interpolation = num( ratio )
    decimation    = den( ratio )
    inLen = ( outputlength * decimation + initialφ - 1 ) / interpolation
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
    inLen = inputlength( outputlength, kernel.ratio, kernel.φIdx )
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
    φStep         = mod( decimation, interpolation )
    φNext         = currentphase + φStep
    φNext         = φNext > interpolation ? φNext - interpolation : φNext
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
    Nφ                 = self.kernel.Nφ
    tapsPerφ           = self.kernel.tapsPerφ
    xLen               = length( x )
    bufLen             = length( buffer )
    reqDlyLineLen      = self.reqDlyLineLen
    outLen             = outputlength( self, xLen )
    criticalYidx       = min( reqDlyLineLen*interpolation, outLen )

    bufLen >= outLen || error( "length( buffer ) must be >= interpolation * length(x)")

    inputIdx = 1
    φ        = 1

    for yIdx in 1:criticalYidx

        accumulator = zero(T)

        for k in 1:tapsPerφ-inputIdx
            @inbounds accumulator += pfb[k, φ] * dlyLine[k+inputIdx-1]
        end

        for k in 1:inputIdx
            @inbounds accumulator += pfb[tapsPerφ-inputIdx+k, φ] * x[k]
        end

        @inbounds buffer[yIdx]  = accumulator
        (φ, inputIdx) = φ == Nφ ? ( 1, inputIdx+1 ) : ( φ+1, inputIdx )
    end

    for yIdx in criticalYidx+1:outLen

        accumulator = zero(T)

        for k in 1:tapsPerφ
            @inbounds accumulator += pfb[ k, φ ] * x[ inputIdx - tapsPerφ + k ]
        end

        @inbounds buffer[yIdx]  = accumulator
        (φ, inputIdx) = φ == Nφ ? ( 1, inputIdx+1 ) : ( φ+1, inputIdx )
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

    outLen = outputlength( xLen-kernel.inputDeficit+1, kernel.ratio, kernel.φIdx )
    bufLen >= outLen || error( "buffer is too small" )

    pfb::PFB{T}        = kernel.pfb
    dlyLine::Vector{T} = self.dlyLine
    interpolation      = num( kernel.ratio )
    decimation         = den( kernel.ratio )
    φIdxStepSize       = mod( decimation, interpolation )
    criticalφIdx       = kernel.Nφ - φIdxStepSize

    inputIdx           = kernel.inputDeficit
    yIdx               = 0

    while inputIdx <= xLen

        accumulator = zero( T )
        yIdx       += 1

        if inputIdx < kernel.tapsPerφ
            hIdx = 1
            for k in inputIdx:self.reqDlyLineLen
                @inbounds accumulator += pfb[ hIdx, kernel.φIdx ] * dlyLine[ k ]
                hIdx += 1
            end

            for k in 1:inputIdx
                @inbounds accumulator += pfb[ hIdx, kernel.φIdx ] * x[ k ]
                hIdx += 1
            end
        else
            hIdx = 1
            for k in inputIdx-kernel.tapsPerφ+1:inputIdx
                @inbounds accumulator += pfb[ hIdx, kernel.φIdx ] * x[ k ]
                hIdx += 1
            end
        end

        buffer[ yIdx ] = accumulator

        inputIdx   += ifloor( ( kernel.φIdx + decimation - 1 ) / interpolation )
        kernel.φIdx = nextphase( kernel.φIdx, kernel.ratio )
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
#       ____ ___ ____ ___ ____ _    ____ ____ ____    ____ _ _    ___          #
#       [__   |  |__|  |  |___ |    |___ [__  [__     |___ | |     |           #
#       ___]  |  |  |  |  |___ |___ |___ ___] ___]    |    | |___  |           #
#==============================================================================#

function filt( h::Vector, x::Vector, ratio::Rational = 1//1 )
    self = FIRFilter( h, ratio )
    filt( self, x )
end
