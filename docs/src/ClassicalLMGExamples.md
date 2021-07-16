# Examples for ClassicalDicke

```@setup examples
push!(LOAD_PATH,"../../src")
on_github=get(ENV, "CI", nothing) == "true"
using Dicke
```

## Exponential growth of the FOTOC using the TWA
In Ref. [Pilatowsky2020](@cite), it was shown that the fidelity out-of-time order 
correlator (FOTOC) corresponding the the unstable fixed point ``(Q,P)=(0,0)`` at the excited-state 
quantum phase transition of the Lipkin-Meshkov-Glick (LMG) model grows eponentially, 
even though it is a regular system. In this example we calculate such quantity 
using the Truncated Wigner Approximation (TWA)

```@example examples
using Dicke.TruncatedWignerApproximation, Dicke.ClassicalLMG
using Plots

fixed_point = ClassicalLMG.Point(Q=0, P=0)

systemLMG = ClassicalLMGSystem(Ω=1, ξ=-1)
j = 500
W = coherent_Wigner_SU2(fixed_point, j =j)
ts=0:0.1:50
if !on_github #hide
  ts=0:10 #hide
end #hide
FOTOC= sum.(
        variance(systemLMG,
        observable = [:Q,:P], 
        N = 5000, 
        ts = ts, 
        distribution = W)
    )

plot(ts,FOTOC,
    yscale=:log10,
    xlabel="time",
    ylabel="FOTOC",
    key=false)
savefig("FOTOCLMG.svg"); nothing #hide
```
![](FOTOCLMG.svg)