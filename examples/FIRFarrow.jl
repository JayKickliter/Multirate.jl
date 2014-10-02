using Multirate
using DSP
using PyPlot

N𝜙              = 32                                           # Number of polyphase partitions
tapsPer𝜙        = 10                                           # N𝜙 * tapsPer𝜙 will be the length of out prototype filter taps
resampleRate    = float64(pi)                                  # Can be any arbitrary resampling rate
polyorder       = 4                                            # Our taps will tranformed into
Nx              = 40                                           # Number of signal samples
t               = 0:Nx-1                                       # Time range
xƒ1             = 0.15                                         # First singal frequency
xƒ2             = 0.3                                          # Second signal frequency
x               = cos(2*pi*xƒ1*t) .+ 0.5sin(2*pi*xƒ2*t*pi)     # Create the two signals and add them
tx              = [0:length(x)-1]                              # Signal time vector
cutoffFreq      = min( 0.45/N𝜙, resampleRate/N𝜙 )              # N𝜙 is also the integer interpolation, so set cutoff frequency accordingly
hLen            = tapsPer𝜙*N𝜙                                  # Total number of filter taps
h               = firdes( hLen, cutoffFreq, DSP.kaiser ) .* N𝜙 # Generate filter taps and scale by polyphase interpolation (N𝜙)
myfilter        = FIRFilter( h, resampleRate, N𝜙, polyorder )  # Construct a FIRFilter{FIRFarrow} object
y               = filt( myfilter, x )                          # Filter x
ty              = [0:length(y)-1]./resampleRate - tapsPer𝜙/2   # Create y time vector. Accout for filter delay so the plots line up

figure(num=nothing, figsize=(6, 6/golden), dpi=100, facecolor="w", edgecolor="k" )
hold(true)
plt.suptitle( "Farrow Filter Resampling, ratio = $(resampleRate)" )
subplot( 211 )
plot( tx, x, "b")
stem( tx, x, linefmt = "b-", markerfmt = "b." )
xlim( 0, maximum(ty) )
xlabel( "Time" )
ylabel( "Amplitude" )
xticks([])
yticks([])
subplot( 212 )
stem( ty, y, linefmt = "r-", markerfmt = "r." )
plot( ty, y, "r-" )
xlim( 0, maximum(ty) )
xlabel( "Time" )
ylabel( "Amplitude" )
xticks([])
yticks([])
hold(false)
savefig("Farrow.png", dpi=100)
