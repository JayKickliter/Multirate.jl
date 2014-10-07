using Multirate
import Radio: wgn
using GtkInteract
using Winston

function powerspectrum( x::Vector, window::Function = blackman )
    xLen = length( x )
    x  .*= window( xLen )
    10*log10(fftshift(abs2(fft( x ))))
end

function channelizerplots( signal::Vector, channelizedSignals::Matrix )
    (samplesPerChannel,Nchannels) = size( channelizedSignals )
    # Create plot tables to hold original spectrum plus spectrum and signal traces of each channel
    table              = Table(2,1)
    subTable           = Table(2,Nchannels)

    # Compute the spectrum of the original signal and add it to the table
    signalSpectrum     = powerspectrum( signal )
    p                  = plot( linspace(-0.5,0.5,length(signalSpectrum)), signalSpectrum )
    setattr( p, "title", "Original Spectrum" )
    table[1,1]         = p
    ylim(-10,60)

    freqs = linspace( -0.5, 0.5, samplesPerChannel )
    t     = linspace( 0, Nchannels*samplesPerChannel-1, samplesPerChannel )
    
    # Plot the spectrum/time-domain traces of each channel and add them to the subtable
    for channel in 1:Nchannels
        thisSignal          = channelizedSignals[:,channel]
        spectrum            = powerspectrum( thisSignal )
        sp                  = plot( freqs, spectrum ); ylim(-10,60)
        setattr( sp, "title", "Channel $channel" )
        tp                  = plot( t, real(thisSignal), "-b", t, imag(thisSignal), "-r" )
        if Nchannels > 3
            setattr( sp.y1, "ticklabels", [])
            setattr( sp.x1, "ticklabels", [])
            setattr( tp.y1, "ticklabels", [])
            setattr( tp.x1, "ticklabels", [])
        end
        subTable[1,channel] = sp
        subTable[2,channel] = tp
    end

    # Position the subtable traces below the main spectrum trace
    table[2,1] = subTable
    
    return table
end



Tx = Complex128 # Datatype for x
ƒs = 1.0        # Input sample rate

@manipulate for ƒ in -0.5:0.05:0.5, Nchannels in [1:9], samplesPerChannel in [128,256,512,1024]

    # Construct a complex sinusoidal signal at frequency ƒ
    signal             = Tx[ complex(cos(2*pi*ƒ*t),sin(2*pi*ƒ*t)) for t in 0:Nchannels*samplesPerChannel-1 ]
    signal            += wgn( length(signal), power=0.1, returnComplex=true )

    # Instantiate a channelizer with Nchannels
    channelizer        = Channelizer( Nchannels, 32 )
    channelizedSignals = filt( channelizer, signal )
   
    # Create the table of plots
    table = channelizerplots( signal, channelizedSignals )

    display(table)
end