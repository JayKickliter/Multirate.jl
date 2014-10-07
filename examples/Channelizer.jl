using Multirate
import Radio: wgn
using GtkInteract
using Winston

function powerspectrum( x::Vector, window::Function = blackman )
    xLen = length( x )
    x  .*= window( xLen )
    10*log10(fftshift(abs2(fft( x ))))
end

Tx = Complex128 # Datatype for x
ƒs = 1.0        # Input sample rate

@manipulate for ƒ in -0.5:0.05:0.5, Nchannels in [2,4,8], samplesPerChannel in [128,256,512,1024]

    # Construct a complex sinusoidal signal at frequency ƒ
    signal             = Tx[ complex(cos(2*pi*ƒ*t),sin(2*pi*ƒ*t)) for t in 0:Nchannels*samplesPerChannel-1 ]
    signal            += wgn( length(signal), 0.1 )

    # Instantiate a channelizer with Nchannels
    channelizer        = Channelizer( Nchannels, 32 )
    channelizedSignals = filt( channelizer, signal )

    # Create plot tables to hold original spectrum plus spectrum and signal traces of each channel
    table              = Table(2,1)
    subTable           = Table(2,Nchannels)

    # Compute the spectrum of the original signal and add it to the table
    signalSpectrum     = powerspectrum( signal )
    p                  = plot( linspace(-0.5,0.5,length(signalSpectrum)), signalSpectrum )
    table[1,1]         = p
    ylim(-10,60)

    # Plot the spectrum/time-domain traces of each channel and add them to the subtable
    for channel in 1:Nchannels
        thisSignal          = channelizedSignals[:,channel]
        t                   = 0:samplesPerChannel-1
        spectrum            = powerspectrum( thisSignal )
        sp                  = plot( spectrum ); ylim(-10,60)
        subTable[1,channel] = sp
        tp                  = plot( t, real(thisSignal), "-b", t, imag(thisSignal), "-r" )
        subTable[2,channel] = tp
    end

    # Position the subtable traces below the main spectrum trace
    table[2,1] = subTable

    display(table)
end
