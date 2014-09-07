reload( "Multirate.jl" )
MR = Multirate
using Base.Test

# http://cnx.org/contents/30153d7b-7130-4dd3-89af-53125f3891e5@10

Nfilters        = 32;
resampleRate    = 0.95
cutoffFreq      = 0.45;
transitionWidth = 0.05;
(hLen, β)       = MR.kaiserlength( transitionWidth, samplerate = Nfilters );
hLen            = iceil(  hLen/Nfilters  ) .* Nfilters;
h               = MR.firdes( hLen, cutoffFreq, MR.kaiser, samplerate = 32, beta = β ) .* Nfilters;


arb_resampler  = MR.FIRFilter( h, resampleRate );
arb_piecewise_resampler = MR.FIRFilter( h, resampleRate );

x = [ 1.0:int(rand(97:110)) ];

@printf( "\n\n\nNaive arbitrary resampling\n" )
@time naivere_result  = MR.naivefilt( h, x, resampleRate, Nfilters);

@printf( "\n\n\nNormal arbitrary resampling\n" )
@time arb_result      = MR.filt( arb_resampler, x );


@printf( "\n\n\nPiecewise arbitrary resampling\n" )
@time begin
    arb_piecewise_result = MR.filt( arb_piecewise_resampler, x[1:10] );
    append!( arb_piecewise_result, MR.filt( arb_piecewise_resampler, x[11:end] ) );
end

display( [ naivere_result arb_result arb_piecewise_result ] )

@test areApprox( naivere_result, arb_result )
@test areApprox( naivere_result, arb_piecewise_result )
