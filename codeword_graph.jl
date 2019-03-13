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


numcwbit = 750
# numsim = Int(1e6)
numsim = 0
numiniter = 1
nummaxiter = 5
iswritefile = false
modprint = 20
jldpath = "./results/tmp.jld"

# G_path = "./H/N_12_R_0.5_InReg_2_G"
# H_path = "./H/N_12_R_0.5_InReg_2_HwithG"
G_path = "./H/N_480_R_0.5_Reg_3_G" ;
H_path = "./H/N_480_R_0.5_Reg_3_HwithG" ;

rng = MersenneTwister(0) #use "randseed:0"
deg = [1, 2]
pmf = [0.346, 0.654 ] # from "A New Design of Raptor Codes for Small Message Length - Yuan"

numcwbit = 14*48*2

C = KRaptor.Raptor(rng, G_path, H_path, deg, pmf, numiniter, nummaxiter)

cw1 = []
cw2 = []

for _ in 1:10
	global cw1,cw2

	KRaptor.setnewLT!(C, numcwbit)
	data = bitrand(rng, 240 )
	cw = KRaptor.encode(C, data)
	cw1 = vcat(cw1,cw)
	cw2 = vcat(cw2,cw[1:15*44*2])
end

println(size(cw2))

writedlm( "cw1.csv",  cw1, ',')
writedlm( "cw2.csv",  cw2, ',')