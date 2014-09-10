resamp     = 1.5
δ_hr       = 0.0
δ_acc      = 0.0
δ_lr       = 0.0
count      = 0
numfilters = 32
∇          = inv(resamp)
tau        = 0.0

while count < 10
    count = count + 1
    tau  += ∇
    hIdx = floor( tau * numfilters )
    count % numfilters == 0 && (tau -= 1)
    println( "$count    $tau, $hIdx")
end
