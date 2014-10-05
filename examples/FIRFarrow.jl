using Multirate
using DSP
using PyPlot

                                                                 # function main()
N𝜙                = 32                                           # Number of polyphase partitions
tapsPer𝜙          = 8                                            # N𝜙 * tapsPer𝜙 will be the length of out protoyTimepe filter taps
ƒsIn              = 1.0                                          # Input sample rate
ƒsOut             = π                                            # Input sample rate
resampleRatio     = ƒsOut/ƒsIn
polyorder         = 4                                            # Our taps will tranformed into
xƒ1               = 0.125                                        # First singal frequency
xƒ2               = 0.3                                          # Second signal frequency
xLen              = 20                                           # Number of signal samples
xTime             = [0:xLen-1]                                   # Signal time vector
x                 = cos(2*pi*xƒ1*xTime)
x                 = x + 0.5sin(2*pi*xƒ2*xTime*pi)                # Create the two signals and add them
x                 = x + cos(0.1*xTime)
cutoffFreq        = min( 0.45/N𝜙, resampleRatio/N𝜙 )             # N𝜙 is also the integer interpolation, so set cutoff frequency accordingly
hLen              = tapsPer𝜙*N𝜙                                  # Tintal number of filter taps
h                 = firdes( hLen, cutoffFreq, DSP.kaiser ) .* N𝜙 # Generate filter taps and scale by polyphase interpolation (N𝜙)
myfilter          = FIRFilter( h, resampleRatio, N𝜙, polyorder ) # Construct a FIRFilter{FIRFarrow} object
kernel            = myfilter.kernel

δout              = 3.5                                          # Specified absolute time offset for first output sample
δfilter           = (hLen-1)/(2*N𝜙)                              # The time delay (at output rate) caused by the filtering process
(phase,throwaway) = modf(δfilter+δout)                           # Number of samples (at output rate) to throw away to meet time offset 𝜙

kernel.inputDeficit += throwaway
setphase( myfilter, phase )

y                 = filt( myfilter, x )                          # Filter x
yTime             = [0.0:length(y)-1]/resampleRatio + δout       # Create y time vector. Accout for filter delay so the plots line up




figure(num=1, figsize=(10, 10/golden), dpi=100, facecolor="w", edgecolor="k" )
clf()
rc( "font", size=10 )
hold(true)
plt.suptitle( "Farrow Filter Resampling, ratio = $(resampleRatio)" )

subplot( 311 )

plot( xTime, x, "b")
stem( xTime, x, linefmt = "b-", markerfmt = "b." )
xlabel( "Time" )
ylabel( "Amplitude" )

subplot( 312 )
stem( yTime, y, linefmt = "r-", markerfmt = "r." )
plot( yTime, y, "r-" )
xlabel( "Time" )
ylabel( "Amplitude" )
xlim( yTime[1], xTime[end]+yTime[1] )

subplot( 313 )
plot( xTime, x, "b.-")
plot( yTime, y, "r.-" )

hold(false)
