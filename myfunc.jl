module MyFunc

using Random

function andsum(x :: BitVector, y :: BitVector)
    if size(x) != size(y)
        throw(ArgumentError("Dimensions don't match"))
    end
    tmp = x .& y
    sum = false
    for elem in tmp
        sum = sum ⊻ elem
    end
    return sum
end

function andsum(x :: BitVector, A :: BitMatrix)
    # used for encoding, like out = xA
    # x: input data vector
    # A: generator matrix
    # out: output codeword

    if length(x) != size(A,1)
        throw(ArgumentError("Dimensions don't match"))
    end

    out = BitVector(undef,size(A,2))
    for n = 1 : size(A,2)
        out[n] = andsum(x, A[:,n])
    end
    return out
end

function BPSK(cw :: BitVector)
    # Mapping:
    #   0 -> +1
    #   1 -> -1
    out = ComplexF64.(cw)
    out = out .* (-2.0) .+ 1.0
    return out
end

function QPSK(cw :: BitVector)

    b1 = BPSK(cw[1:2:end])
    b2 = BPSK(cw[2:2:end])

    out = b1 + b2*im
    return out
end


function BPSK(cw :: Vector{Int})
    # Mapping:
    #   0 -> +1
    #   1 -> -1
    #  -1 -> 0

    numofbit = length(cw)
    out = Vector{ComplexF64}(undef, numofbit)
    for n = 1 : numofbit
        if cw[n] == 0
            out[n] = 1.0 + 0.0im
        elseif cw[n] == 1
            out[n] = -1.0 + 0.0im
        elseif cw[n] == -1
            out[n] = 0.0 + 0.0im
        else
            throw("Error input for BPSK !")
        end
    end
    return out
end

mutable struct AWGNC
    rng :: MersenneTwister
    noise_var_all :: Float64 # includes real and imag
    noise_var :: Float64 # only real or imag
    scaling :: Float64

    function AWGNC(rng :: MersenneTwister, noise_var_all :: Float64)
        noise_var = 0.5 * noise_var_all
        scaling = sqrt(noise_var_all)
        return new(rng, noise_var_all, noise_var, scaling)
    end
end

function (chan :: AWGNC)(x :: Number)
    return x + randn(chan.rng, ComplexF64) * chan.scaling
end

function (chan :: AWGNC)(x :: Array)
    return x .+ randn(chan.rng, ComplexF64, size(x)) .* chan.scaling
end


function demodBPSK(noise_var :: Float64, y :: Number)
    # noise_var: only for real part
    return 2.0 * real(y) / noise_var
end

function demodBPSK(noise_var :: Float64, y :: Array)
    # noise_var: only for real part
    return 2.0 .* real(y) ./ noise_var
end

function demodQPSK(noise_var :: Float64, y :: Array)
    # noise_var: only for real part

    # why *2 ?
    b1_LLR = 2.0 .* real(y) ./ noise_var
    b2_LLR = 2.0 .* imag(y) ./ noise_var
    bits_LLR = hcat(b1_LLR,b2_LLR)'
    dim = prod(size(bits_LLR))
    bits_LLR = reshape(bits_LLR,dim)
end

function demodQPSK_map_reverse(noise_var :: Float64, y :: Array)
    # noise_var: only for real part
    b1_LLR = -2.0 .* real(y) ./ noise_var[1:2:end]
    b2_LLR = -2.0 .* imag(y) ./ noise_var[2:2:end]
    bits_LLR = hcat(b1_LLR,b2_LLR)'
    dim = prod(size(bits_LLR))
    bits_LLR = reshape(bits_LLR,dim)
end


mutable struct ErrorCounter
    infoperblk :: Int # number of information bits per block
    bitcount :: Int # counter for total bits
    errbitcount :: Int # counter for error bits
    blkcount :: Int # counter for total blocks
    errblkcount :: Int # counter for error blocks
    ber :: Float64 # bit error rate
    bler :: Float64 # block error rate

    function ErrorCounter(infoperblk :: Int)
        return new(infoperblk, 0, 0, 0, 0, -1.0, -1.0)
    end

end

function (errcount :: ErrorCounter)(data :: BitVector, rec_data :: BitVector)
    if (length(data) != errcount.infoperblk) || (length(rec_data) != errcount.infoperblk)
        throw("Data length not the same !")
    end

    errbitcount_before = errcount.errbitcount
    for n = 1 : errcount.infoperblk
        if data[n] != rec_data[n]
            errcount.errbitcount += 1
        end
    end

    if errbitcount_before != errcount.errbitcount
        errcount.errblkcount += 1
    end

    errcount.bitcount += errcount.infoperblk
    errcount.blkcount += 1

    errcount.ber = errcount.errbitcount / errcount.bitcount
    errcount.bler = errcount.errblkcount / errcount.blkcount
end

function (errcount :: ErrorCounter)(command :: String)
    if command == "reset!"
        errcount.bitcount = 0
        errcount.errbitcount = 0
        errcount.blkcount = 0
        errcount.errblkcount = 0
        errcount.ber = -1.0
        errcount.bler = -1.0
    else
        throw("Wrong command !")
    end
end

function print(errcount :: ErrorCounter)
    println("blkcount = $(errcount.blkcount)")
    println("errblkcount = $(errcount.errblkcount)")
    println("BLER = $(errcount.ber)")
    println("BLER = $(errcount.bler)")
end

function todB(x :: Float64)
    return 10.0 * log10(x)
end

function invdB(x :: Float64)
    return 10.0 ^ (x/10.0)
end

macro getname(arg)
   string(arg)
end

# mutable struct WriteSet
#     dic :: Dict{String,Any}
#     function WriteSet()
#         dic = Dict{String,Any}()
#         return new(dic)
#     end
# end

# function push!(wrset :: Dict{String, Any}, x)
#     nameofx = @getname(x)
#     wrset[nameofx] = x
# end

function calcconf(prob, numsample)
    # calculate the confidence interval
    # sigma2 is the one-side interval with 95 % confidence
    sigma2 = 2.0 * sqrt( prob * (1.0-prob) / numsample)
    return sigma2
end

function calcsample(prob, sigma2 = 0.1 * prob)
    # calculate the necessary number of samples for a given confidence interval
    numsample = ceil(Int64, (( prob * (1.0-prob) ) / ( (sigma2 * 0.5) ^ 2.0 )) )
    return numsample
end

end  # module MyFunc
