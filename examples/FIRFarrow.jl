using Multirate
using DSP
using PyPlot

N𝜙           = 32                                               # Number of polyphase partitions
tapsPer𝜙     = 10                                               # N𝜙 * tapsPer𝜙 will be the length of out protoyTimepe filter taps
resampleRate = float64(pi)                                      # Can be any arbitrary resampling rate
polyorder    = 4                                                # Our taps will tranformed into
xƒ1          = 0.15                                             # First singal frequency
xƒ2          = 0.3                                              # Second signal frequency
xLen         = 40                                               # Number of signal samples
xTime        = [0:xLen-1]                                       # Signal time vector
x            = cos(2*pi*xƒ1*xTime) .+ 0.5sin(2*pi*xƒ2*xTime*pi) # Create the two signals and add them
cutoffFreq   = min( 0.45/N𝜙, resampleRate/N𝜙 )                  # N𝜙 is also the integer interpolation, so set cutoff frequency accordingly
hLen         = tapsPer𝜙*N𝜙                                      # Total number of filter taps
h            = firdes( hLen, cutoffFreq, DSP.kaiser ) .* N𝜙     # Generate filter taps and scale by polyphase interpolation (N𝜙)
myfilter     = FIRFilter( h, resampleRate, N𝜙, polyorder )      # Construct a FIRFilter{FIRFarrow} object
y            = filt( myfilter, x )                              # Filter x
yDelay       = tapsPer𝜙/2                                       # The time delay caused by the filtering process
yTime        = [0:length(y)-1]./resampleRate                    # Create y time vector. Accout for filter delay so the plots line up

figure(num=nothing, figsize=(6, 6/golden), dpi=100, facecolor="w", edgecolor="k" )
rc("font", size=10)
hold(true)
plt.suptitle( "Farrow Filter Resampling, ratio = $(resampleRate)" )
subplot( 211 )
plot( xTime, x, "b")
stem( xTime, x, linefmt = "b-", markerfmt = "b." )
xlim( 0, maximum(xTime)-yDelay )
xlabel( "Time" )
ylabel( "Amplitude" )
subplot( 212 )
stem( yTime, y, linefmt = "r-", markerfmt = "r." )
plot( yTime, y, "r-" )
xlim( yDelay, maximum(yTime) )
xlabel( "Time" )
ylabel( "Amplitude" )
hold(false)
# savefig( "Farrow.png", dpi=100 )