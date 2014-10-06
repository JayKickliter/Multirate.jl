using Multirate
using PyPlot

N𝜙              = 32
tapsPer𝜙        = 10
resampleRate    = 3.14
xƒ              = 0.1
Nx              = 25
x               = cos(2*pi*xƒ*(0:Nx-1))
# x             = zeros( Nx); x[1] = 1.0
cutoffFreq      = min( 0.45, resampleRate )
transitionWidth = 0.05
hLen            = tapsPer𝜙*N𝜙
h               = firdes( hLen, cutoffFreq, DSP.kaiser, samplerate = N𝜙 ) .* N𝜙
myfilter        = FIRFilter( h, resampleRate, N𝜙 )
δfilter         = (hLen-1)/(2*N𝜙)
y               = filt( myfilter, x )
tx              = [0:length(x)-1]
ty              = [0:length(y)-1]./resampleRate-δfilter


stem( tx, x, "r" )
hold(true)
stem( ty, y )
hold(false)
