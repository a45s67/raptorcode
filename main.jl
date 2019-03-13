include("./kldpc.jl")
include("./myfunc.jl")
include("./kraptor.jl")

# using Revise
import .KLDPC
import .MyFunc
import .KRaptor
# using Pkg
# Pkg.add("DelimitedFiles")
using DelimitedFiles
using Random
# using Statistics
using JLD
using Distributed
using Printf

israndseed = true
randseed = 72
EsN0dB_list = collect(-5.0:1:30.0)
numcwbit = 750
# numsim = Int(1e6)
numsim = 0
numiniter = 10
nummaxiter = 5
iswritefile = false
modprint = 20
jldpath = "./results/tmp.jld"

# G_path = "./H/N_12_R_0.5_InReg_2_G"
# H_path = "./H/N_12_R_0.5_InReg_2_HwithG"
G_path = "./H/N_480_R_0.5_Reg_3_G" ;
H_path = "./H/N_480_R_0.5_Reg_3_HwithG" ;

if israndseed == true
    rng = MersenneTwister(randseed)
else
    rng = MersenneTwister()
end

if iswritefile == false
    fpath = "./results/ttttt.txt"
end

deg = [1, 2]
pmf = [0.346, 0.654 ] # from "A New Design of Raptor Codes for Small Message Length - Yuan"

# dist = Raptor.DegDist(rng, deg, pmf)
# vv = dist(Int(1e4))
# basket = zeros(Int64, length(deg))
# for elem in vv
#     for d = 1 : length(deg)
#         if elem == deg[d]
#             basket[d] += 1
#             break
#         end
#     end
# end

C = KRaptor.Raptor(rng, G_path, H_path, deg, pmf, numiniter, nummaxiter)
KRaptor.setnewLT!(C, numcwbit)
rng = MersenneTwister(randseed)
data1 = bitrand(rng, 10)

rng = MersenneTwister(randseed)
data2 = bitrand(rng, 10)
println(data2,data1)


 # only real or imag part

errcount = MyFunc.ErrorCounter(C.numinfobit)
errcount("reset!")
savdic = Dict{String, Any}()
savdic["errorcount"] = errcount


encode_num = collect(100:50:1000)

for sim = 1 : numsim
    println("========sim",sim,"=================")
    for EsN0dB in EsN0dB_list
        
        EsN0 = MyFunc.invdB(EsN0dB)
        N0 = 1.0./EsN0
        chan = MyFunc.AWGNC(rng, N0)
        noise_var = 0.5 * N0
        @printf("EsN0 = %2.2f \n",EsN0)
        
        data = bitrand(rng, C.numinfobit)

        rec_data::BitVector = []

        for numcwbit in encode_num
            global numbit = numcwbit

            KRaptor.setnewLT!(C, numcwbit)
            cw = KRaptor.encode(C, data)
            tx_symb = MyFunc.BPSK(cw)
            rx_symb = chan(tx_symb)
            rx_LLR = MyFunc.demodBPSK(noise_var, rx_symb)
            (succes , rec_data ) = KRaptor.decode!(C, rx_LLR)

            if succes
                
                #local rec_data
                println(size(rec_data))
                @goto escape_label
            end
        end

        @label escape_label
        errcount(data, rec_data)
        
        @printf("encode_num = %d , rate = %f\n",numbit ,240.0/numbit)
        MyFunc.print(errcount)
        save(jldpath, savdic)
        

    end
end
println()
println("Simulation End !")
MyFunc.print(errcount)
save(jldpath, savdic)
