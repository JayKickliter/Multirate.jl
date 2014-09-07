reload( expanduser("~/.julia/v0.3/Multirate/src/Multirate.jl") )
const MR = Multirate
using Base.Test

# http://cnx.org/contents/30153d7b-7130-4dd3-89af-53125f3891e5@10

x          = [ 1.0:int(rand(97:110)) ];
upfactor   = 32;
downfactor = 37;
resampleRate = 0.876

filter_size    = 32;
arb_cutoff     = 0.45;
arb_transition = 0.05;
arb_h          = MR.firdes( arb_cutoff, arb_transition, samplerate = filter_size ) .* filter_size ;
arb_resampler  = MR.FIRFilter( arb_h, resampleRate );

arb_piecewise_resampler = MR.FIRFilter( arb_h, resampleRate );

# rational_cutoff     = 0.45 / downfactor;
# rational_transition = 0.05 / downfactor;
# rational_h          = MR.firdes( rational_cutoff, rational_transition ) .* filter_size;
# rational_resampler  = MR.FIRFilter( arb_h, upfactor//downfactor );

@time naivere_result  = MR.naivefilt( arb_h, x, resampleRate, filter_size);
# @time rational_result = MR.filt( rational_resampler, x );
@time arb_result      = MR.filt( arb_resampler, x );

@time begin
    arb_piecewise_result = MR.filt( arb_piecewise_resampler, x[1:10] );
    append!( arb_piecewise_result, MR.filt( arb_piecewise_resampler, x[11:end] ) );
end


commonLen = min( length(rational_result), length(naivere_result), length(arb_result), length(arb_piecewise_result) )

# display( [ [1:commonLen] rational_result[1:commonLen] naivere_result[1:commonLen] arb_result[1:commonLen] arb_piecewise_result[1:commonLen] ] )
# println( "Lengths: rational_result $(length(rational_result)), naivere_result $(length(naivere_result)), arb_result $(length(arb_result))")
# @test areApprox( naivere_result, rational_result )
display( [ naivere_result arb_result arb_piecewise_result ] )

@test areApprox( naivere_result, arb_result )
@test areApprox( naivere_result, arb_piecewise_result )
