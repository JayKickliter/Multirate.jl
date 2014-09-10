
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
        ∇            = inv(rate)
        𝜙Accumulator = 0.0
        𝜙IdxVirtual  = 1.0
        𝜙Idx         = 1
        Δ            = 0.0
        new( rate, N𝜙, ∇, 𝜙Accumulator, 𝜙Idx, Δ, 𝜙IdxVirtual )
    end
end

function increment!( self::ArbResamplerState )
        self.𝜙Accumulator += self.∇
        self.𝜙IdxVirtual = self.𝜙Accumulator * self.N𝜙
        self.𝜙Idx        = ifloor( self.𝜙IdxVirtual )
        self.Δ           = mod( self.𝜙Idx, 1 )

        if self.𝜙Idx > self.N𝜙
            self.𝜙Accumulator -= 1
            self.𝜙IdxVirtual  -= self.N𝜙
            self.𝜙Idx          = ifloor( self.𝜙IdxVirtual )
            self.Δ             = self.𝜙IdxVirtual - self.𝜙Idx
        end
end


resamp = 0.9
N𝜙     = 32
yCount = 0
xCount = 0
self   = ArbResamplerState( resamp, 32 )

while xCount < 10
    xCount += 1
    if self.𝜙Idx
    println( "$yCount: 𝜙Idx = $(self.𝜙Idx), Δ = $(self.Δ), 𝜙IdxVirtual = $(self.𝜙IdxVirtual), 𝜙Accumulator = $(self.𝜙Accumulator)")
    increment!( self )
end


while xCount < 10
    xCount += 1
    while 𝜙Idx <= N𝜙
        yCount       += 1

        println( "$yCount: 𝜙Idx = $𝜙Idx, Δ = $Δ, 𝜙IdxVirtual = $𝜙IdxVirtual, 𝜙Accumulator = $𝜙Accumulator")

        𝜙Accumulator += ∇
        𝜙IdxVirtual   = 𝜙Accumulator * N𝜙
        𝜙Idx          = ifloor( 𝜙IdxVirtual ) + 1
        Δ             = mod( 𝜙Accumulator, 1 )

    end

    𝜙Accumulator -= 1
    𝜙IdxVirtual  -= N𝜙
    𝜙Idx          = ifloor( 𝜙IdxVirtual )

end
