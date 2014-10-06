using Multirate
using Winston

powerspectrum( x::Vector ) = 10*log10(fftshift(abs2(fft( x ))))

Nchannels         = 4
Nsignals          = 1
samplesPerChannel = 128
channelizer       = Channelizer( Nchannels )

signalFreqs      = rand( Nsignals )-0.5
signalAmplitudes = rand( Nsignals )
compositeSignal  = zeros( Complex64, Nchannels*samplesPerChannel )

for t in 0:length(compositeSignal)-1
    for signalIdx in 1:Nsignals
        ƒ = signalFreqs[signalIdx]
        A = signalAmplitudes[signalIdx]
        compositeSignal[t+1] += A*exp(im*2*pi*ƒ*t)
    end
end

channelizedSignals = filt( channelizer, compositeSignal )
compositeSpectrum  = powerspectrum( compositeSignal )

mainTable      = Table(2,1)
mainTable[1,1] = plot( linspace(-0.5,0.5,length(compositeSpectrum)), compositeSpectrum )
subTable       = Table(1,Nchannels)

for channel in 1:Nchannels
    spectrum            = powerspectrum( channelizedSignals[:,channel] )
    subTable[1,channel] = plot( linspace(-0.5,0.5,length(spectrum)), spectrum )
end

mainTable[2,1] = subTable

display( mainTable )