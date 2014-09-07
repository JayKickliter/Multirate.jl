reload( "Multirate.jl" )
MR = Multirate
using Base.Test

# http://cnx.org/contents/30153d7b-7130-4dd3-89af-53125f3891e5@10

Nfilters        = 32;
resampleRate    = 9.999;
cutoffFreq      = 0.45;
transitionWidth = 0.05;
(hLen, β)       = MR.kaiserlength( transitionWidth, samplerate = Nfilters );
hLen            = iceil(  hLen/Nfilters  ) .* Nfilters;
h               = MR.firdes( hLen, cutoffFreq, MR.kaiser, samplerate = 32, beta = β ) .* Nfilters;
t               = linspace( 0, 40*2*pi, 100)
x               = sin( t )

arb_resampler  = MR.FIRFilter( h, resampleRate );
arb_piecewise_resampler = MR.FIRFilter( h, resampleRate );

@printf( "\n\n\nNaive arbitrary resampling\n" )
@time naivere_result  = MR.naivefilt( h, x, resampleRate, Nfilters);

@printf( "\n\n\nNormal arbitrary resampling\n" )
@time arbitrary_result      = MR.filt( arb_resampler, x );


@printf( "\n\n\nPiecewise arbitrary resampling\n" )
pivot_point = ifloor( length(x)/4 )
@time begin
    piecewise_result = MR.filt( arb_piecewise_resampler, x[1:pivot_point] );
    append!( piecewise_result, MR.filt( arb_piecewise_resampler, x[pivot_point+1:end] ) );
end

display( [ [1:length(naivere_result)] naivere_result arbitrary_result  piecewise_result naivere_result.-arbitrary_result naivere_result.-piecewise_result ] )

@test areApprox( arbitrary_result,  piecewise_result )
@test areApprox( naivere_result,    arbitrary_result )
@test areApprox( naivere_result,    piecewise_result )