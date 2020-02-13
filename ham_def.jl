
@everywhere tlist = 0:par.t:par.ts

@everywhere aq = destroy(q_basis) ⊗ identityoperator(c_basis)    # qubit annihilation op
@everywhere ac = identityoperator(q_basis) ⊗ destroy(c_basis)    # cavity annihilation op

@everywhere Xq = (aq + dagger(aq)); @everywhere Npq = dagger(aq)*aq;
@everywhere Xc = (ac + dagger(ac)); @everywhere Npc = dagger(ac)*ac;

@everywhere Hq = par.ωq*Npq + (par.ah/2)*(Npq-Npq*Npq)    # qubit Duffing oscillator model
@everywhere Hc = par.ωc*Npc    # cavity - harmonic oscillator
@everywhere Hd = Xc    # cavity drive
@everywhere H0 = Hq + Hc + Hqc    #time independent part of Hamiltonian

@everywhere J = [sqrt(par.kappa_c)*ac, sqrt(par.kappa_q)*aq]
@everywhere Jdagger = dagger.(J)
if coupling_type==1
	println("coupling type = -g*(Xq)*(Xc)")
	@everywhere Hqc = -par.g*(Xq)*(Xc) # transverse qubit-cavity coupling

elseif coupling_type==2
	println("coupling type = -g*(Npq)*(Npc)")
	@everywhere Hqc = -par.g*(Npq)*(Npc) # longitudinal coupling

elseif coupling_type==3
	println("coupling type = -(g/4.0)*(Xq*Xq)*(Xc*Xc)")
	@everywhere Hqc = -(par.g/4.0)*(Xq*Xq)*(Xc*Xc) # longitudinal coupling

elseif coupling_type==4:
	println("coupling type = -g/4*(Xq*Xq)*(Xc*Xc) - g_trans*(Xq)*(Xc)")
	@everywhere Hqc = -(par.g/4.0)*(Xq*Xq)*(Xc*Xc) - g_trans*(Xq)*(Xc) # transverse and longitudinal coupling

elseif coupling_type==5:
	println("coupling type = - g * (Xq.cosm() - 1) * (Xc.cosm() -1 )")
	@everywhere Hqc = - par.g * (Xq.cosm()-1) * par.g_trans*(Xc.cosm()-1) # longitudinal coupling

end

@everywhere function H(t, psi) # time-dependent Hamiltonian
    H = H0 + -1*A*cos(ωd*t + ϕ)*Xc
    return H, J, Jdagger
end
