
type ArbResamplerState
    rate::Float64
    N𝜙::Int
    ∇::Float64
    𝜙Accumulator::Float64
    𝜙Idx::Int
    Δ::Float64
    𝜙IdxVirtual::Float64
    function ArbResamplerState( rate::Real, N𝜙::Integer = 32 )
        rate         = rate
        N𝜙           = N𝜙
        ∇            = 1.0/rate
        𝜙Accumulator = 0.0
        𝜙IdxVirtual  = 0.0
        𝜙Idx         = 0.0
        Δ            = 0.0
        new( rate, N𝜙, ∇, 𝜙Accumulator, 𝜙Idx, Δ, 𝜙IdxVirtual )
    end
end

function increment!( self::ArbResamplerState )
        self.𝜙Accumulator += self.∇

        if self.𝜙Accumulator > 1.0
            self.𝜙Accumulator = mod(self.𝜙Accumulator, 1.0) 
        end
        display(self.𝜙Accumulator)
        self.𝜙IdxVirtual = self.𝜙Accumulator * self.N𝜙
        self.𝜙Idx        = ifloor( self.𝜙IdxVirtual )
        self.Δ           = self.𝜙IdxVirtual - self.𝜙Idx
        
        nothing
end


resamp = 10
N𝜙     = 10
yCount = 0
xCount = 0
self   = ArbResamplerState( resamp, N𝜙 )

while xCount < 30
    xCount += 1
    @printf( "%d: \t𝜙Accumulator = %f\t𝜙IdxVirtual = %f\t𝜙Idx = %f\tΔ = %f\n", xCount, self.𝜙Accumulator, self.𝜙IdxVirtual, self.𝜙Idx, self.Δ)
    increment!( self )
end
#
# resamp       = 0.9
# N𝜙           = 32
# yCount       = 0
# xCount       = 0
# 𝜙Idx         = 0
# Δ            = 0.0
# 𝜙IdxVirtual  = 0.0
# 𝜙Accumulator = 0.0
# ∇ = int(resamp)
#
# while xCount < 10
#     xCount += 1
#     while 𝜙Idx <= N𝜙
#         yCount       += 1
#
#         println( "$yCount: 𝜙Idx = $𝜙Idx, Δ = $Δ, 𝜙IdxVirtual = $𝜙IdxVirtual, 𝜙Accumulator = $𝜙Accumulator")
#
#         𝜙Accumulator += ∇
#         𝜙IdxVirtual   = 𝜙Accumulator * N𝜙
#         𝜙Idx          = ifloor( 𝜙IdxVirtual ) + 1
#         Δ             = mod( 𝜙Accumulator, 1 )
#
#     end
#
#     𝜙Accumulator -= 1
#     𝜙IdxVirtual  -= N𝜙
#     𝜙Idx          = ifloor( 𝜙IdxVirtual )
#
# end
