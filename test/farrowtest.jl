using Multirate

# Filter tap parameters
Th              = Float64
Tx              = Float64
N𝜙              = 32
resamplerate    = 32/33.5
ƒcutoff         = 0.45
transitionwidth = 0.1
(hLen, β)       = kaiserlength( transitionwidth, samplerate = N𝜙 )
hLen            = iceil(  hLen/N𝜙  ) .* N𝜙
h               = Multirate.firdes( hLen, ƒcutoff, DSP.kaiser, samplerate = N𝜙, beta = β ) .* N𝜙
h               = convert( Vector{Th}, h )

# Filter instantiation
polyorder       = 5
x               = rand( Tx, 50  )
arbfilt         = FIRFilter( h, resamplerate, N𝜙 )
farrowfilt      = FIRFilter( h, resamplerate, N𝜙, polyorder )

# Results
yarb            = filt( arbfilt, x)
yfarrow         = filt( farrowfilt, x )
diff_arb_farrow = yarb.-yfarrow
maxError        = maxabs( diff_arb_farrow )
minError        = minabs( diff_arb_farrow )

display( [ [1:length(yarb)] yarb yfarrow diff_arb_farrow] )
println( "Max error: $maxError, Min error: $minError" )