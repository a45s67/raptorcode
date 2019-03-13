include("./kldpc.jl")
include("./myfunc.jl")
include("./kraptor.jl")

# using Revise
import .KLDPC
import .MyFunc
import .KRaptor


cw = BitVector([0,0,1,1])


a = MyFunc.QPSK(cw)
b = MyFunc.demodQPSK(1.0,a)
println(a)
println(b)