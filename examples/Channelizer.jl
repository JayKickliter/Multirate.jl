using Multirate
import Radio: wgn
using GtkInteract
using Winston

function powerspectrum( x::Vector )
    x .*= hanning(length(x))
    10*log10(fftshift(abs2(fft( x ))))
end

Tx = Complex128

@manipulate for ƒ in -0.5:0.05:0.5, Nchannels in [2,4,8], samplesPerChannel in [128,256,512,1024]#, out = :plot

    signal             = Tx[ complex(cos(2*pi*ƒ*t),sin(2*pi*ƒ*t)) for t in 0:Nchannels*samplesPerChannel-1 ]
    signal            += wgn(length(signal), 0.1) 

    channelizer        = Channelizer( Nchannels )
    channelizedSignals = filt( channelizer, signal )

    signalSpectrum     = powerspectrum( signal )
    table              = Table(2,1)
    p                  = plot( linspace(-0.5,0.5,length(signalSpectrum)), signalSpectrum )
    table[1,1]         = p; ylim(0,50)    
    subTable           = Table(1,Nchannels)

    for channel in 1:Nchannels
        spectrum            = powerspectrum( channelizedSignals[:,channel] )
        p                   = plot( spectrum ); ylim(0,50)        
        subTable[1,channel] = p
    end

    table[2,1] = subTable
    display(table)
    
    # push!(out, table)
end
