require( "Multirate.jl" )

MR = Multirate
using Base.Test

# http://cnx.org/contents/30153d7b-7130-4dd3-89af-53125f3891e5@10

Nfilters        = 64;
resampleRate    = 0.99999;
cutoffFreq      = 0.45;
transitionWidth = 0.05;
(hLen, β)       = MR.kaiserlength( transitionWidth, samplerate = Nfilters );
hLen            = iceil(  hLen/Nfilters  ) .* Nfilters;
h               = MR.firdes( hLen, cutoffFreq, MR.kaiser, samplerate = 32, beta = β ) .* Nfilters;
t               = linspace( 0, 1, 200 )
f               = 20
x               = cos( 2*pi*f*t )

@printf( "\n\n\nNaive arbitrary resampling\n" )
@time Ynaive  = MR.naivefilt( h, x, resampleRate, Nfilters);

@printf( "\n\n\nNormal arbitrary resampling\n" )
myfilt  = MR.FIRFilter( h, resampleRate, Nfilters );
@time Ystateless      = MR.filt( myfilt, x );


@printf( "\n\n\nPiecewise arbitrary resampling\n" )
myfilt = MR.FIRFilter( h, resampleRate, Nfilters );
pivot_point      = ifloor( length(x)/4 )
Ypiecwise = eltype(x)[]
sizehint( Ypiecwise, iceil( length(x)*resampleRate ) )
@time for i in 1:length( x )
    thisY = MR.filt( myfilt, x[i:i] )
    append!( Ypiecwise, thisY )
    # println( "Loop $i, y = $thisY")
end


@test areApprox( Ystateless, Ypiecwise )
@test areApprox( Ynaive,     Ystateless )
@test areApprox( Ynaive,     Ypiecwise )

commonLen = min( length(Ynaive), length(Ystateless), length(Ypiecwise) )
display( [ [1:commonLen] Ynaive[1:commonLen] Ystateless[1:commonLen]  Ypiecwise[1:commonLen] Ynaive[1:commonLen].-Ystateless[1:commonLen] Ynaive[1:commonLen].-Ypiecwise[1:commonLen] ] )