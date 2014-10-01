module Multirate


include( "enum.jl" )

import Base: filt, filt!
import Polynomials: Poly, polyval

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
        FIRInterpolator,
        FIRArbitrary,
        FIRDecimator,
        FIRFarrow,
        FIRRational,
        FIRStandard,
        filt!,
        filt,
        tapsforphase!,
        tapsforphase,
        taps2pfb,
        reset!,
        outputlength,
        inputlength
        
include( "NaiveResamplers.jl" )

end # module