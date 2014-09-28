using Winston
using Multirate

Th              = Float64
Tx              = Float64
N𝜙              = 32
resamplerate    = 0.987654321
cutofffreq      = 0.45
transitionwidth = 0.05
(hLen, β)       = kaiserlength( transitionwidth, samplerate = N𝜙 )
hLen            = iceil(  hLen/N𝜙  ) .* N𝜙
h               = Multirate.firdes( hLen, cutofffreq, DSP.kaiser, samplerate = N𝜙, beta = β ) .* N𝜙
h               = convert( Vector{Th}, h )
polyorder       = 6
x             = rand( Tx, 1_000 )
# x               = ones( Tx, 1_000_000)
arbfilt         = FIRFilter( h, resamplerate, N𝜙 )
farrowfilt      = FIRFilter( h, resamplerate, N𝜙, polyorder )

filt( h, x[1:1], resamplerate, N𝜙 )
filt( h, x[1:1], resamplerate, N𝜙, polyorder )

@time yarb              = filt( arbfilt, x)
@time (yfarrow,db_𝜙vec) = filt( farrowfilt, x )

diff_arb_farrow = yarb.-yfarrow
maxError        = maxabs( diff_arb_farrow )
minError        = minabs( diff_arb_farrow )

display( [[1:length(yarb)] yarb yfarrow diff_arb_farrow db_𝜙vec] )



println( "Max error: $maxError, Min error: $minError" )
