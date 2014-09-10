resamp     = .85
δ_hr       = 0.0
δ_acc      = 0.0
δ_lr       = 0.0
count      = 0
numfilters = 32
Δ          = 1/numfilters

while count < 100
    count = count + 1

    δ_hr  = mod( (count-1) * numfilters / resamp, 1 )
    δ_lr = δ_lr + Δ
    δ_lr  = δ_lr >= 1.0 ? δ_lr-1.0 : δ_lr
    
    println( "$count: δ_hr = $δ_hr, δ_lr = $δ_lr")
end
