using QuantumOptics

#=
    Input:
        fq : frequency of Qubit
        fc : frequency of Cavity
        fd:  cavity drive frequency
        ah1: Anharmonicity of Qubit
        g_qc : Qubit-Cavity coupling
        kappa_c: line width of the cavity

        Cavity drive of type  A1*cos(wt+phi)
        A1: Amplitude of drive
        phi: phase of drive
        t:time length of drive
        ts:number of point per time unit.
	    t_cutoff:time after which drive stops

    Output:
        defines following global variables like,
        aq,Npq,ac,Npc : annihilation and number operator of Q1
        H0     : wq*Nq + (ah/2)*Nq*(1-Nq) + w1*Nq1
        H1     : (a + a.dag())
        c_op_list: list containing collapse operators iff Tq is non-zero
                   Only 1-->0 decays is considered. (Sinlge photon decay)

        returns density matrix for given tlist
=#

function instate(Nq,q_state_amp,Nc)
#Nq = levels in the qubit
#state_amp = amplitudes of the state
#Nc = levels in the cavity
    if mod(Nq,2) == 0
        q_basis = SpinBasis((Nq-1) // 2)
    else
        q_basis = SpinBasis((Nq-1)/2 // 1)
    end

    c_basis = FockBasis(Nc)
    ψ0 = Ket(q_basis, q_state_amp) ⊗ fockstate(c_basis,0)

    return q_basis, c_basis, ψ0
end

function H(t, ρ) # time-dependent Hamiltonian
    if t<t_cutoff
        H = H0 + A*cos(ωd*t + ϕ)*Iqc
    else
        H = H0
    end

    return H, J, Jdagger
end

q_basis, c_basis, ψ0 = instate(6,[0,1,0,0,0,0],20)

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
ωd = 2.0*pi*6.78     # drive frequency
ϕ = 0                # phase of drive
t_cutoff = 60       # drive cutoff time

aq = sigmam(q_basis) ⊗ identityoperator(c_basis)    # qubit annihilation op
ac = identityoperator(q_basis) ⊗ destroy(c_basis)    # cavity annihilation op

Iqc = identityoperator(q_basis) ⊗ identityoperator(c_basis) # Identity operator

Xq = (aq + dagger(aq)); Yq = -1*im*(aq - dagger(aq)); Zq = (Iqc - 2*dagger(aq)*aq); Npq = dagger(aq)*aq;
Xc = (ac + dagger(ac)); Yc = -1*im*(ac - dagger(ac)); Zc = (Iqc - 2*dagger(ac)*ac); Npc = dagger(ac)*ac;

Hq = ωq*Npq + (ah/2)*Npq*(Iqc-Npq)    # qubit Duffing oscillator model
Hc = ωc*Npc    # cavity - harmonic oscillator
Hqc = g*(Xq)*(Xc)    # qubit cavity interaction
Hd = Xc    # cavity drive
H0 = Hq + Hc + Hqc    #time independent part of Hamiltonian

J = [sqrt(kappa_c)*ac,sqrt(kappa_q)*aq]
Jdagger = dagger.(J)
#R = [kappa_c,kappa_q]
T = 1:0.025:100
@time tout, ρ = timeevolution.master_dynamic(T, ψ0, H)
