using Distributed
if(nprocs()==1)
    addprocs(2)
end

@everywhere using QuantumOptics, DifferentialEquations, JLD2, SharedArrays

@everywhere function instate(Nq,q_state_amp,Nc)
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

@everywhere function H(t, psi) # time-dependent Hamiltonian
    if t<t_cutoff
        H = H0 + -1*A*cos(ωd*t + ϕ)*Xc
    else
        H = H0
    end

    return H, J, Jdagger
end

@everywhere nq = 5
@everywhere nc = 9
@everywhere q_basis, c_basis, ψ0 = instate(nq,[0,1,0,0,0,0],nc)

@everywhere fq = 5.4
@everywhere fc = 6.78
@everywhere fd = 6.78
@everywhere ah1 = 0.221
@everywhere g_qc = 0.087
@everywhere kappa_c = 0.027
@everywhere kappa_q = 0.02

@everywhere ωq = fq*2.0*pi       # qubit freq
@everywhere ωc = fc*2.0*pi       # cavity freq
@everywhere ah = ah1*2.0*pi      # qubit anharmonicity
@everywhere g  = g_qc*2.0*pi     # qubit cavity coupling const
@everywhere A  = 2.0*pi*0.1      # drive amplitude
@everywhere ωd = 2.0*pi*fd     # drive frequency
@everywhere ϕ = 0                # phase of drive
@everywhere t_cutoff = 60       # drive cutoff time

@everywhere aq = destroy(q_basis) ⊗ identityoperator(c_basis)    # qubit annihilation op
@everywhere ac = identityoperator(q_basis) ⊗ destroy(c_basis)    # cavity annihilation op

@everywhere Xq = (aq + dagger(aq)); @everywhere Npq = dagger(aq)*aq;
@everywhere Xc = (ac + dagger(ac)); @everywhere Npc = dagger(ac)*ac;

@everywhere Hq = ωq*Npq + (ah/2)*(Npq-Npq*Npq)    # qubit Duffing oscillator model
@everywhere Hc = ωc*Npc    # cavity - harmonic oscillator
@everywhere Hqc = g*(Xq)*(Xc)    # qubit cavity interaction
@everywhere Hd = Xc    # cavity drive
@everywhere H0 = Hq + Hc + Hqc    #time independent part of Hamiltonian

@everywhere J = [sqrt(kappa_c)*ac, sqrt(kappa_q)*aq]
@everywhere Jdagger = dagger.(J)
@everywhere T = 0:0.01:100

println(Threads.nthreads())
#@time tout, ρ = timeevolution.master_dynamic(T, ψ0, H)
#@save "data "*string(nq)*","*string(nc)*","*string(maximum(T))*".jld2" tout ρ T

@everywhere Ntrajectories = 10

@time @sync @distributed for i = 1:Ntrajectories
    t_tmp, ψ = timeevolution.mcwf_dynamic(T, ψ0, H; seed=i)
    id = myid()
    println(id)
    f = open("trajectories_from_$id.txt","a")
    write(f, "$ψ \n -- \n")
    close(f)
    println(i,"th seed complete")
end

#@save "data "*string(nq)*","*string(nc)*","*string(maximum(T))*".jld2" Ntrajectories T

rmprocs(workers())
println("Successfully removed workers")
