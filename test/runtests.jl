using Base.Test
import Multirate
import DSP

# Disable time and printf macros when not running interactivly ( for travis )
if !isinteractive()
    macro time( ex )
        quote
            $(esc(ex))
        end
    end

    macro printf( args... )
    end
end

function areApprox( x1::Vector, x2::Vector )
    Nx1 = length( x1 )
    Nx2 = length( x2 )

    if Nx1 != Nx2
        # @printf("x1 & x2 are different lengths vectors")
        return false
    end

    for i = 1:Nx1
        if !isapprox( x1[i], x2[i] )
            # @printf( "Something went wrong at index $i" )
            return false
        end
    end

    return true
end

function naiveresample{T}( h::Vector{T}, x::Vector{T}, resampleRate::Real, numfilters::Integer = 32 )

    xLen = length( x )
    yLen = iceil( xLen * resampleRate )
    y    = similar( x, yLen )
    x    = filt( h, x, numfilters//1 )

    yIdx = 1

    for yIdx in 1:yLen
        xIdxVirtual = numfilters * (yIdx - 1) / resampleRate + 1
        xIdxLower   = ifloor( xIdxVirtual )
        xIdxUpper   = iceil( xIdxVirtual )
        Δ           = xIdxVirtual - xIdxLower
        y[yIdx]     = x[xIdxLower] + Δ*( x[xIdxUpper] - x[xIdxLower] )
    end

    return y
end

#==============================================================================#
#               ____ _ _  _ ____ _    ____    ____ ____ ___ ____               #
#               [__  | |\ | | __ |    |___    |__/ |__|  |  |___               #
#               ___] | | \| |__] |___ |___    |  \ |  |  |  |___               #
#==============================================================================#

function test_singlerate( h, x )
    xLen       = length( x )
    hLen       = length( h )
    pivotPoint = min( rand(50:150, 1)[1], ifloor( xLen/4 ))
    x1         = x[ 1 : pivotPoint ]
    x2         = x[ pivotPoint+1 : end ]

    @printf( "\n\n" )
    @printf( "____ _ _  _ ____ _    ____    ____ ____ ___ ____\n" )
    @printf( "[__  | |\\ | | __ |    |___    |__/ |__|  |  |___\n" )
    @printf( "___] | | \\| |__] |___ |___    |  \\ |  |  |  |___\n" )
    @printf( "\nTesting single-rate fitering, h is %s, x is %s. xLen = %d, hLen = %d", string(eltype(h)), string(eltype(x)), xLen, hLen )

    @printf( "\n\tBase.filt\n\t\t")
    @time baseResult = Base.filt( h, 1.0, x )

    if method_exists( DSP.firfilt, ( typeof(h), typeof(x) ))
        @printf( "\n\tDSP.firfilt\n\t\t")
        @time dspResult = DSP.firfilt( h, x )
    end

    @printf( "\n\tMultirate.filt( h, x, 1//1 )\n\t\t" )
    @time statelesResult = Multirate.filt( h, x )

    @printf( "\n\tMultirate.filt. length( x1 ) = %d, length( x2 ) = %d\n\t\t", length( x1 ), length( x2 ) )
    self = Multirate.FIRFilter( h, 1//1 )
    @time begin
        y1 = Multirate.filt( self, x1 )
        y2 = Multirate.filt( self, x2 )
    end
    statefulResult = [ y1, y2 ]

    @printf( "\n\tMultirate.filt filt. Piecewise for first %d inputs\n\t\t", length( x1 ) )
    Multirate.reset( self )
    @time begin
        for i in 1:length(x1)
            y1[i] = Multirate.filt( self, x1[i:i] )[1]
        end
        y2 = Multirate.filt( self, x2 )
    end
    piecewiseResult = [ y1, y2 ]


     if areApprox( baseResult, statelesResult ) && areApprox( baseResult, statefulResult ) && areApprox( baseResult, piecewiseResult )
         return true
     end

     display( [ baseResult statelesResult statefulResult piecewiseResult ] )

     return false

end




#==============================================================================#
#                      ___  ____ ____ _ _  _ ____ ___ ____                     #
#                      |  \ |___ |    | |\/| |__|  |  |___                     #
#                      |__/ |___ |___ | |  | |  |  |  |___                     #
#==============================================================================#

function test_decimation( h, x, decimation )
    xLen       = length( x )
    hLen       = length( h )
    pivotPoint = min( rand(50:150, 1)[1], ifloor( xLen/4 ))
    x1         = x[ 1 : pivotPoint ]
    x2         = x[ pivotPoint+1 : end ]

    @printf( "\n\n" )
    @printf( "___  ____ ____ _ _  _ ____ ___ _ ____ _  _ \n" )
    @printf( "|  \\ |___ |    | |\\/| |__|  |  | |  | |\\ | \n" )
    @printf( "|__/ |___ |___ | |  | |  |  |  | |__| | \\| \n" )
    @printf( "\nTesting decimation. h::%s, x::%s. xLen = %d, hLen = %d, decimation = %d", string(typeof(h)), string(typeof(h)), xLen, hLen, decimation )

    @printf( "\n\tNaive decimation with Base.filt\n\t\t")
    @time begin
        baseResult   = Base.filt( h, one(eltype(h)), x )
        baseResult   = baseResult[1:decimation:end]
    end

    if method_exists( DSP.firfilt, ( typeof(h), typeof(x) ))
        @printf( "\n\tNaive decimation with DSP.firfilt\n\t\t")
        @time begin
            dspResult = DSP.firfilt( h, x )
            dspResult = dspResult[1:decimation:end]
        end
    end

    @printf( "\n\tMultirate.filt( h, x, 1//%d)\n\t\t", decimation )
    @time statelesResult = Multirate.filt( h, x, 1//decimation )

    @printf( "\n\tMultirate.filt decimation. length( x1 ) = %d, length( x2 ) = %d\n\t\t", length( x1 ), length( x2 ) )
    self = Multirate.FIRFilter( h, 1//decimation )
    @time begin
        y1 = Multirate.filt( self, x1 )
        y2 = Multirate.filt( self, x2 )
    end
    statefulResult = [ y1, y2 ]

    @printf( "\n\tMultirate.filt decimation. Piecewise for first %d inputs.\n\t\t", length( x1 ) )
    Multirate.reset( self )
    y1 = similar( x, 0 )
    @time begin
        for i in 1:length(x1)
            append!( y1, Multirate.filt( self, x1[i:i] ) )
        end
        y2 = Multirate.filt( self, x2 )
    end
    piecewiseResult = [ y1, y2 ]

    if areApprox( baseResult, statelesResult ) && areApprox( baseResult, statefulResult ) && areApprox( baseResult, piecewiseResult )
        return true
    end

    display( [ baseResult statefulResult piecewiseResult ] )
    return false

end




#==============================================================================#
#               _ _  _ ___ ____ ____ ___  _    ____ ____ ___ ____              #
#               | |\ |  |  |___ |__/ |__] |    |  | |__|  |  |___              #
#               | | \|  |  |___ |  \ |    |___ |__| |  |  |  |___              #
#==============================================================================#

function test_interpolation( h, x, interpolation )
    xLen       = length( x )
    hLen       = length( h )
    pivotPoint = min( rand(50:150, 1)[1], ifloor( xLen/4 ))
    x1         = x[ 1 : pivotPoint ]
    x2         = x[ pivotPoint+1 : end ]

    @printf( "\n\n" )
    @printf( "_ _  _ ___ ____ ____ ___  _    ____ ____ ___ _ ____ _  _ \n" )
    @printf( "| |\\ |  |  |___ |__/ |__] |    |  | |__|  |  | |  | |\\ | \n" )
    @printf( "| | \\|  |  |___ |  \\ |    |___ |__| |  |  |  | |__| | \\| \n" )
    @printf( "\nTesting interpolation, h::%s, x::%s. xLen = %d, hLen = %d, interpolation = %d", typeof(h), typeof(x), xLen, hLen, interpolation )

    @printf( "\n\tNaive interpolation with Base.filt\n\t\t")
    @time begin
        xZeroStuffed = zeros( eltype(x), xLen * interpolation )
        for n = 0:xLen-1;
            xZeroStuffed[ n*interpolation+1 ] = x[ n+1 ]
        end
        baseResult = Base.filt( h, one(eltype(h)), xZeroStuffed )
    end

    if method_exists( DSP.firfilt, ( typeof(h), typeof(x) ))
        @printf( "\n\tNaive interpolation with DSP.firfilt\n\t\t")
        @time begin
            xZeroStuffed = zeros( eltype(x), xLen * interpolation )
            for n = 0:xLen-1;
                xZeroStuffed[ n*interpolation+1 ] = x[ n+1 ]
            end
            dspResult = DSP.firfilt( h, xZeroStuffed )
        end
    end

    @printf( "\n\tMultirate.filt( h, x, %d//1 )\n\t\t", interpolation )
    @time statelesResult = Multirate.filt( h, x, interpolation//1 )

    @printf( "\n\tMultirate.filt interpolation. length( x1 ) = %d, length( x2 ) = %d\n\t\t", length( x1 ), length( x2 ) )
    self = Multirate.FIRFilter( h, interpolation//1 )
    @time begin
        y1 = Multirate.filt( self, x1 )
        y2 = Multirate.filt( self, x2 )
    end
    statefulResult = [ y1, y2 ]

    @printf( "\n\tMultirate.filt interpolation. Piecewise for first %d inputs\n\t\t", length( x1 ) )
    Multirate.reset( self )
    y1 = similar( x, 0 )
    @time begin
        for i in 1:length(x1)
            append!( y1, Multirate.filt( self, x1[i:i] ) )
        end
        y2 = Multirate.filt( self, x2 )
    end
    piecewiseResult = [ y1, y2 ]

    if areApprox( baseResult, statelesResult ) && areApprox( baseResult, statefulResult ) && areApprox( piecewiseResult, baseResult )
        return true
    end

    display( [ baseResult statefulResult piecewiseResult ] )

end




#==============================================================================#
#           ____ ____ ___     ____ ____ ____ ____ _  _ ___  _    ____          #
#           |__/ |__|  |      |__/ |___ [__  |__| |\/| |__] |    |___          #
#           |  \ |  |  |  .   |  \ |___ ___] |  | |  | |    |___ |___          #
#==============================================================================#

function test_rational( h, x, ratio )
    xLen       = length( x )
    hLen       = length( h )
    pivotPoint = min( rand(50:150, 1)[1], ifloor( xLen/4 ))
    x1         = x[ 1 : pivotPoint ]
    x2         = x[ pivotPoint+1 : end ]
    upfactor   = num( ratio )
    downfactor = den( ratio )
    resultType = promote_type( eltype(h), eltype(x) )

    @printf( "\n\n" )
    @printf( "      ____ ____ ___ _ ____ _  _ ____ _    \n" )
    @printf( "      |__/ |__|  |  | |  | |\\ | |__| |    \n" )
    @printf( "      |  \\ |  |  |  | |__| | \\| |  | |___ \n" )
    @printf( "                                          \n" )
    @printf( "____ ____ ____ ____ _  _ ___  _    _ _  _ ____\n" )
    @printf( "|__/ |___ [__  |__| |\\/| |__] |    | |\\ | | __\n" )
    @printf( "|  \\ |___ ___] |  | |  | |    |___ | | \\| |__]\n" )
    @printf( "\n\nTesting rational resampling, h::%s, x::%s. xLen = %d, hLen = %d, ratio = %d//%d", string(typeof(h)), string(typeof(x)), xLen, hLen, upfactor, downfactor )

    @printf( "\n\tNaive rational resampling with Base.filt\n\t\t")
    @time begin
        xStuffed   = zeros( resultType, length(x) * upfactor )
        baseResult = Array( resultType, int( ceil( length(x) * ratio )))

        for n = 0:length(x)-1;
            xStuffed[ n*upfactor+1 ] = x[ n+1 ]
        end

        baseResult = Base.filt( h, one(eltype(h)), xStuffed );
        baseResult = [ baseResult[n] for n = 1:downfactor:length( baseResult ) ]
    end

    if method_exists( DSP.firfilt, ( typeof(h), typeof(x) ))
        @printf( "\n\tNaive rational resampling DSP.firfilt\n\t\t")
        @time begin
            xStuffed  = zeros( resultType, length(x) * upfactor )
            dspResult = Array( resultType, int( ceil( length(x) * ratio )))

            for n = 0:length(x)-1;
                xStuffed[ n*upfactor+1 ] = x[ n+1 ]
            end

            dspResult = DSP.firfilt( h, xStuffed );
            dspResult = [ dspResult[n] for n = 1:downfactor:length( dspResult ) ]
        end
    end

    @printf( "\n\tMultirate.filt( h, x, %d//%d )\n\t\t", upfactor, downfactor )
    @time statelesResult = Multirate.filt( h, x, ratio )

    @printf( "\n\tMultirate.filt rational resampling. length( x1 ) = %d, length( x2 ) = %d\n\t\t", length( x1 ), length( x2 ) )
    self = Multirate.FIRFilter( h, ratio )
    @time begin
        s1 = Multirate.filt( self, x1 )
        s2 = Multirate.filt( self, x2 )
    end
    statefulResult = [ s1, s2 ]

    @printf( "\n\tMultirate.filt rational. Piecewise for all %d inputs\n\t\t", length( x ) )
    self = Multirate.FIRFilter( h, ratio )
    y1 = similar( x, 0 )
    @time begin
        for i in 1:length(x)
            append!( y1, Multirate.filt( self, x[i:i] ) )
        end
    end
    piecewiseResult = y1

    if areApprox( baseResult, statelesResult ) && areApprox( baseResult, statefulResult ) && areApprox( baseResult, piecewiseResult )
        return true
    end

    display( [  baseResult statefulResult statefulResult piecewiseResult ] )

    return false
end




function test_all()
    for interpolation in 1:16,
            decimation in 1:16,
                Th in [Float32, Float64],
                    Tx in [Float32, Float64, Complex64, Complex128]

        h     = rand(Th, rand(16:128,1)[1] )
        xLen  = int(rand( 1000:2000, 1 )[1])
        xLen  = xLen-mod( xLen, decimation )
        x     = rand( Tx, xLen )
        ratio = interpolation//decimation

        if ratio == 1
            @test test_singlerate( h, x )
        end

        if decimation != 1
            @test test_decimation( h, x, decimation )
        end

        if interpolation != 1
            @test test_interpolation( h, x, interpolation )
        end


        if num(ratio) != 1 && den(ratio) != 1
            @test test_rational( h, x, ratio )
        end
    end
end

function test_nextphase()
    for interpolation in 1:8
        for decimation in 1:8
            ratio           = interpolation//decimation
            interpolation   = num(ratio)
            decimation      = den(ratio)
            x               = repmat( [1:interpolation], decimation )
            reference = [ x[n] for n = 1:decimation:length( x ) ]
            result = [ 1 ]
            for i in 2:interpolation
                append!( result, [ Multirate.nextphase( result[end], ratio ) ] )
            end
            @test areApprox( reference, result )
        end
    end
end

test_nextphase()
test_all()
