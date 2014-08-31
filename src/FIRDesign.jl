include( "enum.jl" )
include( "Window.jl" )

#==============================================================================#
#               ____ ____ _  _ ____ ___ ____ _  _ ___ ____                     #
#               |    |  | |\ | [__   |  |__| |\ |  |  [__                      #
#               |___ |__| | \| ___]  |  |  | | \|  |  ___]                     #                      
#==============================================================================#

@enum( FIRResponse, LOWPASS, BANDPASS, HIGHPASS, BANDSTOP )




#==============================================================================#
#                                Kaiser Order                                  #
#==============================================================================#
# δ  = ripple in passband and stopband
# Δω = ω_p - ω_s
#    = transition width in radians 
# 
# returns:
#  Ntaps = length of filter
#        
function kaiserord( stopbandAttenuation::Real, transitionWidth::Real, sampleRate = 1.0 )
    
    A     = stopbandAttenuation
    Δω    = 2π * ( transitionWidth/sampleRate )
    Ntaps = iceil( ( A - 7.95 )/( 2.285 * Δω )) + 1
    
    if A > 50
        β = 0.1102*( A - 8.7 )
    elseif A > 21
        β = 0.5842*( A - 21 )^( 0.4 ) + 0.07886*( A - 21 )
    else
        β     = 0.0
        Ntaps = iceil(  5.79 / Δω ) + 1
    end

    return Ntaps, β
end




#==============================================================================#
#          ____ _ ____    ___  ____ ____ ___ ____ ___ _   _ ___  ____          #
#          |___ | |__/    |__] |__/ |  |  |  |  |  |   \_/  |__] |___          #
#          |    | |  \    |    |  \ |__|  |  |__|  |    |   |    |___          #
#==============================================================================#
# Ntaps       = Desired number of filter taps
#             = If FIRType is HIGHPASS, may return  samples to make it
#               a type 1 filter.
# F           = Cutoff frequency. Real for high-pass & lowpass, Vector{Real} for
#               band-pass & band-reject
# FIRResponse = The response of the filter: LOWPASS, BANDPASS, HIGHPASS, BANDSTOP

function firprototype( Ntaps::Integer, F::Union(Real, Vector), response::FIRResponse = LOWPASS )
    M = Ntaps-1
    if response == LOWPASS
        prototype = [ 2*F*sinc(2*F*(n-M/2)) for n = 0:M ]
    elseif response == BANDPASS
        prototype = [ 2*(F[1]*sinc(2*F[1]*(n-M/2)) - F[2]*sinc(2*F[2]*(n-M/2))) for n = 0:M ]        
    elseif response == HIGHPASS
        M = isodd( M ) ? M+1 : M 
        prototype = [ sinc(n-M/2) - 2*F*sinc(2*F*(n-M/2)) for n = 0:M ]        
    elseif response == BANDSTOP
        prototype = [ 2*(F[2]*sinc(2*F[2]*(n-M/2)) - F[1]*sinc(2*F[1]*(n-M/2))) for n = 0:M ]
    else
        error("Not a valid FIR_TYPE")
    end
    
    return prototype
end




#==============================================================================#
#                        ____ _ ____ ___  ____ ____                            #
#                        |___ | |__/ |  \ |___ [__                             #
#                        |    | |  \ |__/ |___ ___]                            #
#==============================================================================#

function firdes( Ntaps::Integer, cutoff::Union(Real, Vector), windowFunction::Function; response::FIRResponse = LOWPASS, sampleRate = 1.0, beta = 5.0 )
    
    cutoff    = cutoff ./ sampleRate
    prototype = firprototype( Ntaps, cutoff, response )

    if windowFunction == kaiser              
        return prototype .* kaiser( Ntaps, beta ) 
    else
        return prototype .*  windowFunction( Ntaps )
    end

end
