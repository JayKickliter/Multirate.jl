# Borrowed from http://www.mathworks.com/help/signal/ref/upfirdn.html

using  Winston
import Multirate

sampleRate    = 48000
interpolation = 147
decimation    = 160
ratio         = interpolation  // decimation
Nt            = 24*interpolation

h = Multirate.firdes( Nt, 0.5/interpolation, Multirate.kaiser, beta = 7.8562  )

n = 0:10239
x = sin(2*pi*1e3/sampleRate*n)
y = Multirate.filt( h, x, ratio )

tOriginal = n[1:49]./sampleRate
xOriginal = x[1:49]
tResamp   = n[1:45]./(sampleRate*ratio)
xResamp   = y[13:57]

display( stem( tOriginal, xOriginal, "b." ) )
hold( true )
display( stem( tResamp, xResamp, "r." ) )
hold( false )