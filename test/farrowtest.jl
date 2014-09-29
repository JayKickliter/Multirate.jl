using Winston
using Multirate

function plotphase( farrowfilt::FIRFilter, phase::Real )
    farrowkernel = farrowfilt.kernel
    pntaps       = [ Polynomials.polyval( farrowkernel.pnfb[tapIdx], phase  ) for tapIdx in 1:farrowkernel.tapsPer𝜙 ]
    p            = plot( pntaps )
    display( p )
end

function plotphase( farrowfilt::FIRFilter, arbfilt::FIRFilter, phase::Real, showdelta = false )

    𝜙Idx         = ifloor( phase )
    α            = phase-𝜙Idx
    farrowkernel = farrowfilt.kernel
    arbkernel    = arbfilt.kernel
    tapsPer𝜙     = farrowkernel.tapsPer𝜙
    farrowtaps   = [ Polynomials.polyval( farrowkernel.pnfb[tapIdx], phase  ) for tapIdx in 1:tapsPer𝜙 ]
    arbtaps      = arbkernel.pfb[:,𝜙Idx] .+ α*arbkernel.dpfb[:,𝜙Idx]

    if showdelta
        tapsdelta   = arbtaps.-farrowtaps
        p           = plot( tapsdelta )
        display( p )
    else
        x = 1:tapsPer𝜙
        p = plot( x, arbtaps, x, farrowtaps )
        display( p )
    end
end

function plotrow( farrowfilt::FIRFilter, tapIdx::Real, showdelta = false )

    farrowkernel = farrowfilt.kernel
    pntaps       = [ Polynomials.polyval( farrowkernel.pnfb[tapIdx], phase  ) for phase in 1:farrowkernel.N𝜙 ]
    if showdelta
        x        = 1:farrowkernel.N𝜙
        basetaps = vec(farrowkernel.pfb[tapIdx,:])
        delta    = basetaps.-pntaps
        p        = plot( x, delta )
    else
        p = plot( pntaps )
        display( p )
    end
end


Th              = Float64
Tx              = Float64
N𝜙              = 32
resamplerate    = 32/33.5
cutofffreq      = 0.45
transitionwidth = 0.05
(hLen, β)       = kaiserlength( transitionwidth, samplerate = N𝜙 )
hLen            = iceil(  hLen/N𝜙  ) .* N𝜙
h               = Multirate.firdes( hLen, cutofffreq, DSP.kaiser, samplerate = N𝜙, beta = β ) .* N𝜙
h               = convert( Vector{Th}, h )
polyorder       = 6
x             = rand( Tx, 50  )
# x               = ones( Tx, 1_000_000)
arbfilt         = FIRFilter( h, resamplerate, N𝜙 )
farrowfilt      = FIRFilter( h, resamplerate, N𝜙, polyorder )

filt( h, x[1:1], resamplerate, N𝜙 )
filt( h, x[1:1], resamplerate, N𝜙, polyorder )

(yarb, dbvxarb, dbphiarb) = filt( arbfilt, x)
(yfarrow, dbvxf, dbphif)  = filt( farrowfilt, x )

diff_arb_farrow = yarb.-yfarrow
maxError        = maxabs( diff_arb_farrow )
minError        = minabs( diff_arb_farrow )

display( [ [1:length(yarb)] yarb yfarrow diff_arb_farrow  dbvxarb dbvxf dbphiarb dbphif] )
println( "Max error: $maxError, Min error: $minError" )
