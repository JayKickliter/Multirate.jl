using Multirate
using DSP
using PyPlot

N𝜙              = 32
tapsPer𝜙        = 10
resampleRate    = 3.14
xƒ1             = 0.15
xƒ2             = 0.3
Nx              = 40
t               = 0:Nx-1
x               = ( cos(2*pi*xƒ1*t) .+ 0.5sin(2*pi*xƒ2*t*pi) )  #.* hamming( Nx )
tx              = [0:length(x)-1]
cutoffFreq      = min( 0.45, resampleRate )
transitionWidth = 0.05
hLen            = tapsPer𝜙*N𝜙
h               = firdes( hLen, cutoffFreq, DSP.kaiser, samplerate = N𝜙, beta = 5 ) .* N𝜙
myfilter        = FIRFilter( h, resampleRate, N𝜙, 4 )
y               = filt( myfilter, x )
ty              = [0:length(y)-1]./resampleRate - tapsPer𝜙/2

clf()
hold(true)
p1 = subplot( 211 )
# title( "Original signal")
plot( tx, x, "b-")
xlim( 0, maximum(ty) )
stem( tx, x, "b-o" )
subplot( 212 )
# title( "Resampled at $(resampleRate)x with a Farrow filter")
stem( ty, y, "r-o" )
plot( ty, y, "r-" )
xlim( 0, maximum(ty) )
hold(false)
savefig("Farrow.png", dpi=100)