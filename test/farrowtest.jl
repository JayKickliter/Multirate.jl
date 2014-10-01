using Multirate

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
x               = rand( Tx, 50  )
arbfilt         = FIRFilter( h, resamplerate, N𝜙 )
farrowfilt      = FIRFilter( h, resamplerate, N𝜙, polyorder )
yarb            = filt( arbfilt, x)
yfarrow         = filt( farrowfilt, x )
diff_arb_farrow = yarb.-yfarrow
maxError        = maxabs( diff_arb_farrow )
minError        = minabs( diff_arb_farrow )

display( [ [1:length(yarb)] yarb yfarrow diff_arb_farrow] )
println( "Max error: $maxError, Min error: $minError" )