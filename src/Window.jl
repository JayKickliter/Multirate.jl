#==============================================================================#
#                                Blackman Window                               #
#==============================================================================#

function blackman( n::Integer, numTaps::Integer )
    # TODO: add argument valdiation
    α = 0.42
    β = 0.5
    α - β*cos(2*π*n/(numTaps-1)) + (β-α)*cos(4*π*n/(numTaps-1))
end

function blackman( numTaps::Integer )
    # TODO: add argument valdiation
    [ blackman(n, numTaps) for n = 0:(numTaps-1) ]
end

#==============================================================================#
#                                Hamming Window                                #
#==============================================================================#

function hamming( n::Integer, numTaps::Integer )
    # TODO: add argument valdiation
    α = 0.54
    β = 0.46 # 1 - α
    α - β*cos(2*π*n/(numTaps-1))
end

function hamming( numTaps::Integer )
    # TODO: add argument valdiation
    [ hamming(n, numTaps) for n = 0:(numTaps-1) ]
end

#==============================================================================# 
#                                  Hann Window                                 #
#==============================================================================#

function hann( n::Integer, numTaps::Integer )
    # TODO: add argument valdiation
    α = 0.5
    β = 0.5
    α - β*cos(2*π*n/(numTaps-1))
end

function hann( numTaps::Integer )
    # TODO: add argument valdiation
    [ hann(n, numTaps) for n = 0:(numTaps-1) ]
end

#==============================================================================#
#                                Kaiser Window                                 #
#==============================================================================#

function kaiser( n::Integer, numTaps::Integer, β::Real )
    # TODO: add argument valdiation
    num = besseli(0, β*sqrt(1-(2*n/(numTaps-1)-1).^2))
    dem = besseli(0, β)
    num/dem
end

function kaiser( numTaps::Integer, β::Real )
    # TODO: add argument valdiation
    [ kaiser( n, numTaps, β ) for n = 0:(numTaps-1) ]
end