
type ArbResamplerState
    rate::Float64
    N𝜙::Int
    Δ::Float64
    counter::Float64
    𝜙Accumulator::Float64
    𝜙Idx::Int
    δ::Float64
    𝜙IdxVirtual::Float64
    yLower::Number
    function ArbResamplerState( rate::Real, N𝜙::Integer = 32 )
        rate         = rate
        N𝜙           = N𝜙
        Δ            = 1.0/rate
        counter      = 0
        𝜙Accumulator = 0
        𝜙IdxVirtual  = 1
        𝜙Idx         = 1
        δ            = 0
        yLower       = 0
        self = new( rate, N𝜙, Δ, counter, 𝜙Accumulator, 𝜙Idx, δ, 𝜙IdxVirtual, yLower )
        # update( self )
        self
    end
end

function update( self::ArbResamplerState )
    self.counter += self.Δ
    self.𝜙Accumulator += self.Δ
    if self.𝜙Accumulator >= 1
        self.𝜙Accumulator = mod( self.𝜙Accumulator, 1 )
        self.counter += 1
    end

    self.𝜙IdxVirtual = self.𝜙Accumulator * ( self.N𝜙 - 1 ) + 1
    self.𝜙Idx        = ifloor( self.𝜙IdxVirtual )
    self.δ           = self.𝜙IdxVirtual - self.𝜙Idx
    self
end


resamp = 2
N𝜙     = 10
yCount = 0
xCount = 0
self   = ArbResamplerState( resamp, N𝜙 )

while xCount < 11
    xCount += 1
    @printf( "%d:\tcounter = %f, 𝜙Accumulator = %f, 𝜙IdxVirtual = %f, 𝜙Idx = %d, δ = %f\n", xCount, self.counter, self.𝜙Accumulator, self.𝜙IdxVirtual, self.𝜙Idx, self.δ)
    update( self )
end
#
# resamp       = 0.9
# N𝜙           = 32
# yCount       = 0
# xCount       = 0
# 𝜙Idx         = 0
# δ            = 0.0
# 𝜙IdxVirtual  = 0.0
# 𝜙Accumulator = 0.0
# Δ = int(resamp)
#
# while xCount < 10
#     xCount += 1
#     while 𝜙Idx <= N𝜙
#         yCount       += 1
#
#         println( "$yCount: 𝜙Idx = $𝜙Idx, δ = $δ, 𝜙IdxVirtual = $𝜙IdxVirtual, 𝜙Accumulator = $𝜙Accumulator")
#
#         𝜙Accumulator += Δ
#         𝜙IdxVirtual   = 𝜙Accumulator * N𝜙
#         𝜙Idx          = ifloor( 𝜙IdxVirtual ) + 1
#         δ             = mod( 𝜙Accumulator, 1 )
#
#     end
#
#     𝜙Accumulator -= 1
#     𝜙IdxVirtual  -= N𝜙
#     𝜙Idx          = ifloor( 𝜙IdxVirtual )
#
# end
