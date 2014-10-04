using Multirate
using DSP
using PyPlot

function time_firfarrow{Th,Tx}( self::FIRFilter{FIRFarrow{Th}}, x::Vector{Tx} )
    kernel = self.kernel
    xLen   = length( x )
    println( "\nFIRFarrow speed test" )
    @printf( "\tresampling rate  %f\n", resampleRate )
    # @printf( "\tpolynomial order %d\n", polyorder )
    @printf( "\tx type           %s\n", string(Tx) )
    @printf( "\tx length         %d\n", xLen )
    @printf( "\th type           %s\n", string(Th) )
    @printf( "\th length         %d\n", kernel.tapsPer𝜙*kernel.N𝜙 )
    @printf( "\tN𝜙               %d\n", kernel.N𝜙 )
    @printf( "\ttaps per 𝜙       %d\n", kernel.tapsPer𝜙 )
    ( y, elapsed, allocated, z ) = @timed filt( self, x )
    @printf( "\telapsed time (s) %1.3f\n", elapsed )
    @printf( "\tinput samples/s  %1.3e\n", xLen/elapsed )
    @printf( "\toutput samples/s %1.3e\n", length(y)/elapsed )
end


function time_firarbitrary{Th,Tx}( self::FIRFilter{FIRArbitrary{Th}}, x::Vector{Tx} )
    println( "\nFIRArbitrary Speed Test" )
    @printf( "\tresampling rate  %f\n", resampleRate )
    @printf( "\tx type           %s\n", string(Tx) )
    @printf( "\tx length         %d\n", xLen )
    @printf( "\th type           %s\n", string(Th) )
    @printf( "\th length         %d\n", hLen )
    @printf( "\tN𝜙               %d\n", N𝜙 )
    @printf( "\ttaps per 𝜙       %d\n", tapsPer𝜙 )
    ( y, elapsed, allocated, z ) = @timed filt( self, x )    
    @printf( "\telapsed time (s) %1.3f\n", elapsed )
    @printf( "\tinput samples/s  %1.3e\n", xLen/elapsed )
    @printf( "\toutput samples/s %1.3e\n", length(y)/elapsed )
end

N𝜙           = 32                                               # Number of polyphase partitions
tapsPer𝜙     = 10                                               # N𝜙 * tapsPer𝜙 will be the length of out protoyTimepe filter taps
resampleRate = 1.0                                              # Can be any arbitrary resampling rate
polyorder    = 4                                                # Our taps will tranformed into
Th           = Float32
cutoffFreq   = min( 0.45/N𝜙, resampleRate/N𝜙 )                  # N𝜙 is also the integer interpolation, so set cutoff frequency accordingly
hLen         = tapsPer𝜙*N𝜙                                      # Total number of filter taps
h            = firdes( hLen, cutoffFreq, DSP.kaiser ) .* N𝜙     # Generate filter taps and scale by polyphase interpolation (N𝜙)
xLen         = 10_000_000                                       # Number of signal samples
Tx           = Complex64
x            = rand( Tx, xLen )


farrowfilt = FIRFilter( h, resampleRate, N𝜙, polyorder )      # Construct a FIRFilter{FIRFarrow} object
arbfilt    = FIRFilter( h, resampleRate, N𝜙 )
time_firfarrow( farrowfilt, x )
time_firarbitrary( arbfilt, x )

resampleRate = 1/2.123456789
farrowfilt = FIRFilter( h, resampleRate, N𝜙, polyorder )      # Construct a FIRFilter{FIRFarrow} object
arbfilt    = FIRFilter( h, resampleRate, N𝜙 )
time_firfarrow( farrowfilt, x )
time_firarbitrary( arbfilt, x )