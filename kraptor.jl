module KRaptor

include("./kldpc.jl")

import .KLDPC
using Random

@enum LTChkMode LTCHK_INNER = 1
@enum MsgPassMode MSGPASS_LT = 1 MSGPASS_LDPC = 2


mutable struct DegDist
    # Degree Distribution
    rng :: MersenneTwister
    deg :: Vector{Int} # degree
    pmf :: Vector{Float64} # probability mass function
    pmfsize :: Int # size of deg and pmf
    cdf :: Vector{Float64} # Cumulative Distribution Function

    function DegDist(rng :: MersenneTwister, deg :: Vector{Int}, pmf :: Vector{Float64})
        pmfsize = length(pmf)
        if length(deg) != pmfsize
            throw("Sizes of deg and pmf are different !")
        end

        cdf = Vector{Float64}(undef, pmfsize)
        cdf[1] = pmf[1]
        for n = 2 : pmfsize
            cdf[n] = cdf[n-1] + pmf[n]
        end

        if isapprox(cdf[end], 1.0, rtol = eps()) == false
            throw("Sum of pmf is not 1.0 !")
        end

        return new(rng, deg, pmf, pmfsize, cdf)
    end
end

function (degdist :: DegDist)()
    # generate a sample from degree distribution
    sample = rand(degdist.rng, Float64)
    for n = 1 : degdist.pmfsize
        if sample < degdist.cdf[n]
            return degdist.deg[n]
        end
    end
end

function (degdist :: DegDist)(numrand :: Int64)
    # generate a sample vector from degree distribution
    out = Vector{Int64}(undef, numrand)
    for n = 1 : numrand
        out[n] = degdist()
    end
    return out
end

mutable struct LTChkNode <: KLDPC.Node
    idx :: Int
    degofeql :: Int
    edgetoeqlvec :: Vector{KLDPC.BinEdgeLDPC}
    degofout :: Int
    edgetooutvec :: Vector{KLDPC.Edge}

    function LTChkNode(idx :: Int)
        return new(idx, 0, Vector{KLDPC.BinEdgeLDPC}[], 0, Vector{KLDPC.Edge}[])
    end
end

function (chk :: LTChkNode)(mode :: LTChkMode)
    if mode == LTCHK_INNER # msg passing toward LDPC EqlNodes
        if chk.degofeql == 0 # it means the deg from degdist is zero
            # then nothing needs to be done
            return
        else
            LLRin = Vector{Float64}(undef, chk.degofeql + chk.degofout)
            for n = 1 : chk.degofout
                LLRin[n] = chk.edgetooutvec[n].LLRtochk
            end
            for n = 1 : chk.degofeql
                LLRin[n+chk.degofout] = chk.edgetoeqlvec[n].LLR2chk
            end

            LLRout = KLDPC.boxplusvec(LLRin)

            for n = 1 : chk.degofeql
                chk.edgetoeqlvec[n].LLR2eql = LLRout[n+chk.degofout]
            end
        end
    else
        throw("Wrong mode for LTChkNode msg passing !")
    end
end

mutable struct LTHalfBinEdge <: KLDPC.Edge
    # used for receiving LT codeword from channel
    idx :: Int
    LLRtochk :: Float64
    idxofchk :: Int
    LTHalfBinEdge() = new(0, 0.0, 0)
    LTHalfBinEdge(idx :: Int) = new(idx, 0.0, 0)
end

function cnct!(edge :: KLDPC.BinEdgeLDPC, chk :: LTChkNode, eql :: KLDPC.EqlNode)
    edge.idx_chk = chk.idx
    edge.idx_eql = eql.idx
    push!(chk.edgetoeqlvec, edge) # point to the same edge object
    chk.degofeql += 1
    push!(eql.edgevec2out, edge) # point to the same edge object
    eql.deg_out += 1
end

function cnct!(edge :: LTHalfBinEdge, chk :: LTChkNode)
    edge.idxofchk = chk.idx
    push!(chk.edgetooutvec, edge) #  point to the same edge object
    chk.degofout += 1
end

mutable struct LTGenMat
    # Sparse LT generator matrix
    G :: Vector{Vector{Int}}
    numcwbit :: Int # number of output codeword bit
    numinbit :: Int # number of input bit

    function LTGenMat()
        # empty initialize
        G = []
        numcwbit = 0
        numinbit = 0
        return new(G, numcwbit, numinbit)
    end

    function LTGenMat(G :: Vector{Vector{Int}}, numinbit :: Int)
        # initialize by G
        numcwbit = size(G, 1)
        return new(G, numcwbit, numinbit)
    end
end

function encode(gen :: LTGenMat, inputbit :: BitVector)
    # encode inputbit by LTGenMat
    # output a Vector{Int64}
    # if deg == 0, then output -1

    if gen.numinbit != length(inputbit)
        throw("Different dimension between LTGenMat and inputbit")
    end

    out = zeros(Int64, gen.numcwbit)
    for n = 1 : gen.numcwbit
        if length(gen.G[n]) == 0
            out[n] = -1
        else
            for idx in gen.G[n]
                out[n] = out[n] ⊻ inputbit[idx]
            end
        end
    end

    return out
end


mutable struct Raptor
    rng :: MersenneTwister
    ldpc :: KLDPC.LDPC
    degdist :: DegDist
    numiniter :: Int # inner iteration of each stage
    nummaxiter :: Int
    numinfobit :: Int
    numcwbit :: Int
    rate :: Float64
    LTgen :: LTGenMat
    LTchkvec :: Vector{LTChkNode}
    LTedgevec :: Vector{KLDPC.BinEdgeLDPC}
    numLTedge :: Int
    hedgevec :: Vector{LTHalfBinEdge}

    function Raptor(rng :: MersenneTwister, G_path :: String, H_path :: String, deg :: Vector{Int}, pmf :: Vector{Float64}, numiniter :: Int = 7, nummaxiter :: Int = 10000)
        ldpc = KLDPC.LDPC(G_path, H_path)
        numinfobit = ldpc.num_info_bit
        degdist = DegDist(rng, deg, pmf)
        LTgen = LTGenMat() # empty
        numcwbit = 0
        rate = 0.0
        LTchkvec = Vector{LTChkNode}(undef, 0)
        LTedgevec = Vector{KLDPC.BinEdgeLDPC}(undef, 0)
        numLTedge = 0
        hedgevec = Vector{LTHalfBinEdge}(undef, 0)

        return new(rng, ldpc, degdist, numiniter, nummaxiter, numinfobit, numcwbit, rate, LTgen, LTchkvec, LTedgevec, numLTedge, hedgevec)
    end
end

function setnewLT!(C :: Raptor, numcwbit :: Int)
    # set a new LT codebook for the Raptor code

    # remove the original LT codebook first
    KLDPC.clearoutedge!(C.ldpc)
    hedgevec = Vector{LTHalfBinEdge}(undef, 0)
    # some remove actions are as below, point to new objects

    # start a new LT codebook
    C.numcwbit = numcwbit
    C.rate = C.numinfobit / numcwbit
    C.LTchkvec = Vector{LTChkNode}(undef, numcwbit)
    degvec = C.degdist(numcwbit)
    C.numLTedge = sum(degvec)
    C.LTedgevec = Vector{KLDPC.BinEdgeLDPC}(undef, C.numLTedge)

    for g = 1 : C.numLTedge
        C.LTedgevec[g] = KLDPC.BinEdgeLDPC(g)
    end

    G = Vector{Vector{Int64}}(undef, 0)
    for n = 1 : numcwbit
        push!(G, Int64[])
    end

    idxedge = 1
    for n = 1 : numcwbit
        C.LTchkvec[n] = LTChkNode(n)
        idxvec = unifpick(C.rng, C.ldpc.num_code_bit, degvec[n])
        for idxnode in idxvec
            push!(G[n], idxnode)
            cnct!(C.LTedgevec[idxedge], C.LTchkvec[n], C.ldpc.eqlvec[idxnode])
            idxedge += 1
        end
    end
   # println(G)
    C.LTgen = LTGenMat(G, C.ldpc.num_code_bit)
    return G
end

function encode(C :: Raptor, data :: BitVector)
    # encode data by Raptor Code

    ldpccw = KLDPC.encode(C.ldpc, data)
    out = encode(C.LTgen, ldpccw)
    return out
end

# function init!(C :: Raptor)
#     # used for receiving Raptor codeword from p2p channel
#     C.hedgevec = Vector{LTHalfBinEdge}(undef, C.numcwbit)
#     for n = 1 : C.numcwbit
#         C.hedgevec[n] = LTHalfBinEdge(n)
#         cnct!(C.hedgevec[n], C.LTchkvec[n])
#     end
# end

function rxinit!(C :: Raptor, rx_LLR :: Vector{Float64})
    # used for receiving Raptor codeword from p2p channel
    # clear some edge LLR

    C.hedgevec = Vector{LTHalfBinEdge}(undef, C.numcwbit)

    for n = 1 : C.numcwbit
        C.hedgevec[n] = LTHalfBinEdge(n)
        cnct!(C.hedgevec[n], C.LTchkvec[n])
        C.hedgevec[n].LLRtochk = rx_LLR[n]
    end
    for edge in C.LTedgevec
        edge.LLR2chk = 0.0
        edge.LLR2eql = 0.0
    end
    KLDPC.rxinit!(C.ldpc)
end

function msgpass!(C :: Raptor, mode :: MsgPassMode)
    if mode == MSGPASS_LT
        for eql in C.ldpc.eqlvec
            eql(KLDPC.EQLINIT_LLRIN)
        end

        for iter = 1 : C.numiniter
            for chk in C.LTchkvec
                chk(LTCHK_INNER)
            end
            for eql in C.ldpc.eqlvec
                eql(KLDPC.EQLOUTER)
            end
        end
        for eql in C.ldpc.eqlvec
            eql(KLDPC.EQLINNER_WO_CALC)
        end
    elseif mode == MSGPASS_LDPC
        (iscw, rec_data) = KLDPC.msgpass!(C.ldpc, C.numiniter)
        return (iscw, rec_data)
    else
        throw("Wrong mode for Raptor msgpass!")
    end
end

function decode!(C :: Raptor, rx_LLR :: Vector{Float64})
    # used for receiving Raptor codeword from p2p channel

    rxinit!(C, rx_LLR)
    for iter = 1 : C.nummaxiter
        msgpass!(C, MSGPASS_LT)
        (iscw, rec_data) = msgpass!(C, MSGPASS_LDPC)
        if iscw == true
            return true , rec_data
        end
    end

    # if the code goes through nummaxiter, it means the msg passing doesn't converge, probably the decoded data is not correct.
    # still, the function output the final results

    rec_cw = KLDPC.getcw(C.ldpc)
    return false , KLDPC.getdata(C.ldpc, rec_cw)
end


function unifpick_slow(rng :: MersenneTwister, N :: Int, K :: Int)
    # uniformly pick K Int from 1:N
    # never pick the same Int
    permvec = randperm(rng, N)
    return permvec[1:K]
end

function unifpick(rng :: MersenneTwister, N :: Int, K :: Int)
    # uniformly pick K elements from 1:N
    # never pick the same Int
    candidate = Array(1:N)
    out = Vector{Int64}(undef, K)
    for k = 1 : K
        samp = rand(rng, 1:N-k+1)
        out[k] = candidate[samp]
        deleteat!(candidate, samp)
    end
    return out
end

end
