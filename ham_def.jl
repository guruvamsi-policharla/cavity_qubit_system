
tlist = 0:par.ts:par.t

aq = destroy(q_basis) ⊗ identityoperator(c_basis)    # qubit annihilation op
ac = identityoperator(q_basis) ⊗ destroy(c_basis)    # cavity annihilation op

Xq = (aq + dagger(aq));   Npq = dagger(aq)*aq;
Xc = (ac + dagger(ac));   Npc = dagger(ac)*ac;

Hq = par.ωq*Npq + (par.ah/2)*(Npq-Npq*Npq)    # qubit Duffing oscillator model
Hc = par.ωc*Npc    # cavity - harmonic oscillator
Hd = Xc    # cavity drive

if coupling_type==1
	#println("coupling type = -g*(Xq)*(Xc)")
	Hqc = -par.g*(Xq)*(Xc) # transverse qubit-cavity coupling
elseif coupling_type==2
	#println("coupling type = -g*(Npq)*(Npc)")
	Hqc = -par.g*(Npq)*(Npc) # longitudinal coupling
elseif coupling_type==3
	#println("coupling type = -(g/4.0)*(Xq*Xq)*(Xc*Xc)")
	Hqc = -(par.g/4.0)*(Xq*Xq)*(Xc*Xc) # longitudinal coupling
elseif coupling_type==4
	#println("coupling type = -g/4*(Xq*Xq)*(Xc*Xc) - g_trans*(Xq)*(Xc)")
	Hqc = -(par.g/4.0)*(Xq*Xq)*(Xc*Xc) - g_trans*(Xq)*(Xc) # transverse and longitudinal coupling
elseif coupling_type==5
	#println("coupling type = - g * (Xq.cosm() - 1) * (Xc.cosm() -1 )")
	Hqc = - par.g * (Xq.cosm()-1) * par.g_trans*(Xc.cosm()-1) # longitudinal coupling
end

H0 = Hq + Hc + Hqc    #time independent part of Hamiltonian

J = [sqrt(par.κc)*ac, sqrt(par.κq)*aq]
Jdagger = dagger.(J)

function H(t, psi) # time-dependent Hamiltonian
	H = H0 + -1*par.A*cos(par.ωd*t + par.ϕ)*Xc
	return H, J, Jdagger
end
