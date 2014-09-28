# using Winston
# using Polynomials
using Multirate

N𝜙              = 32
resampleRate    = 1.0
cutoffFreq      = 0.45
transitionWidth = 0.05
(hLen, β)       = kaiserlength( transitionWidth, samplerate = N𝜙 )
hLen            = iceil(  hLen/N𝜙  ) .* N𝜙
h               = Multirate.firdes( hLen, cutoffFreq, DSP.kaiser, samplerate = N𝜙, beta = β )
polyorder       = 4
x               = ones(100)
farrowfilt      = FIRFilter( h, resampleRate, N𝜙, polyorder )
arbfilt         = FIRFilter( h, resampleRate, N𝜙 )

@time yfarrow = filt( farrowfilt, x )
@time yarb    = filt( arbfilt, x)
display( typeof(yfarrow) )
display( typeof(yarb) )
display( [[1:length(yarb)] yarb yfarrow yarb.-yfarrow] )
