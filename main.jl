using QuantumOptics
using DifferentialEquations
using JLD2
using PyPlot

function instate(Nq,q_state_amp,Nc)
#Nq = levels in the qubit
#state_amp = amplitudes of the state
#Nc = levels in the cavity
    #=
    if mod(Nq,2) == 0
        q_basis = SpinBasis((Nq-1) // 2)
    else
        q_basis = SpinBasis((Nq-1)/2 // 1)
    end
    =#
    q_basis = FockBasis(Nq)
    c_basis = FockBasis(Nc)
    ψ0 = Ket(q_basis, q_state_amp) ⊗ fockstate(c_basis,0)

    return q_basis, c_basis, ψ0
end

function H(t, ρ=0) # time-dependent Hamiltonian
    if t<t_cutoff
        H = H0 + -1*A*cos(ωd*t + ϕ)*Xc
    else
        H = H0
    end

    return H, J, Jdagger
end

nq = 5
nc = 9
q_basis, c_basis, ψ0 = instate(nq,[0,1,0,0,0,0],nc)

fq = 5.4
fc = 6.78
fd = 6.78
ah1 = 0.221
g_qc = 0.087
kappa_c = 0.027
kappa_q = 0.02

ωq = fq*2.0*pi       # qubit freq
ωc = fc*2.0*pi       # cavity freq
ah = ah1*2.0*pi      # qubit anharmonicity
g  = g_qc*2.0*pi     # qubit cavity coupling const
A  = 2.0*pi*0.1      # drive amplitude
ωd = 2.0*pi*fd     # drive frequency
ϕ = 0                # phase of drive
t_cutoff = 60       # drive cutoff time

aq = destroy(q_basis) ⊗ identityoperator(c_basis)    # qubit annihilation op
ac = identityoperator(q_basis) ⊗ destroy(c_basis)    # cavity annihilation op

Xq = (aq + dagger(aq)); Npq = dagger(aq)*aq;
Xc = (ac + dagger(ac)); Npc = dagger(ac)*ac;

Hq = ωq*Npq + (ah/2)*(Npq-Npq*Npq)    # qubit Duffing oscillator model
Hc = ωc*Npc    # cavity - harmonic oscillator
Hqc = g*(Xq)*(Xc)    # qubit cavity interaction
Hd = Xc    # cavity drive
H0 = Hq + Hc + Hqc    #time independent part of Hamiltonian

J = [sqrt(kappa_c)*ac, sqrt(kappa_q)*aq]
Jdagger = dagger.(J)
#R = [kappa_c,kappa_q]
T = 0:0.01:100

println(Threads.nthreads())
@time tout, ρ = timeevolution.master_dynamic(T, ψ0, H;alg=OrdinaryDiffEq.Vern7())
@save "data "*string(nq)*","*string(nc)*","*string(maximum(T))*".jld2" tout ρ T
pygui(true)
plot(tout, (expect(Npc, ρ)))
grid("on")
