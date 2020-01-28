using QuantumOptics

#=
main(fq,fc,fd,ah1,g_qc,kappa_c,A1,phi,t,ts,t_cutoff):
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

function main(fq,fc,fd,ah1,g_qc,kappa_c,A1,phi,t,ts,t_cutoff)
    wq = fq*2.0*pi       # qubit freq
    wc = fc*2.0*pi       # cavity freq
    wd = fd*2.0*pi       # drive freq
    ah = ah1*2.0*pi      #qubit anharmonicity
    g  = g_qc*2.0*pi     # qubit cavity coupling const
    A  = A1*2.0*pi      # drive amplitude


end
