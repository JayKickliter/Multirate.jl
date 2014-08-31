# Borrowed from http://www.mathworks.com/help/signal/ref/upfirdn.html

import Multirate
import Winston

sampleRate    = 48000
interpolation = 147
decimation    = 160
ratio         = interpolation  // decimation
Nt            = 24*interpolation

h = Multirate.firdes( Nt, 0.5/decimation, Multirate.kaiser, beta = 7.8562  ) .* interpolation

n = 0:10239                    
x = sin(2*pi*1e3/sampleRate*n) 
y = Multirate.filt( h, x, ratio )

tOriginal = n[1:49]./sampleRate
xOriginal = x[1:49]
tResamp   = n[1:45]./(sampleRate*ratio)
xResamp   = y[13:57]

plot(tOriginal, xOriginal, "b:", tResamp, xResamp, "r:")