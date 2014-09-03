module Multirate

include( "enum.jl" )

import Base: filt, filt!, reset

using DSP.Windows
export  hanning,
        hammming,
        kaiser,
        blackman

include( "FIRDesign.jl" )
export  firdes,
        kaiserord,
        FIRResponse,
        LOWPASS,
        HIGPASS,
        BANDPASS,
        BANDSTOP

include( "Filters.jl" )
export  FIRFilter,
        polyize,
        filt!,
        filt,
        reset,
        outputlength,
        inputlength

end # module