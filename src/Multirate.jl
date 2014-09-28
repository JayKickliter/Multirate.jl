module Multirate

include( "enum.jl" )

import Base: filt, filt!, reset

using DSP.Windows
export  hanning,
        hammming,
        kaiser,
        blackman

include( "support.jl" )
include( "FIRDesign.jl" )
export  firdes,
        kaiserlength,
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

include( "NaiveResamplers.jl" )

end # module