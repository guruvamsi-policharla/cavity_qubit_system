function instate(nq,qub_amp,nc,cav_amp,st_type)
	#Nq = levels in the qubit
	#qub_amp = amplitude vector of qubit state
	#Nc = levels in the cavity
	#st_type = type of state to be prepared. Fock-Fock by default
	#cav_amp = amplitude vector of cavity state/a scalar alpha if coherence state
    q_basis = FockBasis(nq)
    c_basis = FockBasis(nc)
    if (st_type=="FF")
        ψ0 = Ket(q_basis, qub_amp) ⊗ fockstate(c_basis,cav_amp)
    elseif (st_type == "FC")
        ψ0 = Ket(q_basis, qub_amp) ⊗ coherentstate(c_basis,cav_amp)
	else
		error("Basis not found")
    end
    return q_basis, c_basis, ψ0
end

@with_kw mutable struct param_list
            ωq::Float64 = 0
            ωc::Float64 = 0
            ωd::Float64 = 0
            ah::Float64 = 0
            g::Float64 = 0
			g_trans::Float64 = 0
            A::Float64 = 0
            κq::Float64 = 0
            κc::Float64 = 0
            ϕ::Float64 = 0
            t::Float64 = 0
            ts::Float64 = 0
            t_cutoff::Float64 = 0
            qub_amp = 0
            cav_amp = 0
            nq::Int64 = 0
            nc::Int64 = 0
        end

function define_params_transverse()
    par = param_list()
    par.nq = 5
    par.nc = 6
    par.qub_amp=[0,1,0,0,0,0]    #size needs to be nq+1
    par.t = 1000
    par.ts = 0.8
    par.t_cutoff = 800

    par.ωq = 4.47 * 2*pi       # qubit freq
	par.ωc = 7.415 * 2*pi       # cavity freq
	par.ωd = 7.4211 * 2*pi       # drive freq FOR QB IN EXCITED STATE
	#par.ωd = 7.4218*2.0*np.pi       # drive freq FOR QB IN GROUND STATE
	par.ah = 0.211 * 2*pi      #qubit anharmonicity
	par.g = 0.163 * 2*pi     # qubit cavity coupling const
	par.A = 0.0053 * 2*pi      # drive amplitude
	#par.A  = 0.0078 *2*pi      # drive amplitude for about 25 photons in cavity
	if kappa_on
		par.κc = 0.0014 * 2*pi #cavity linewidth
	else
		par.κc = 0.0 * 2*pi #cavity linewidth
	end
	par.κq = 0 *2*pi #qubit linewidth
	par.ϕ = 0
	return par
end

function define_params_longitudinal()
    par = param_list()
    par.nq = 5
    par.nc = 6
    par.qub_amp=[0,1,0,0,0,0]    #size needs to be nq+1
    par.t = 1000
    par.ts = 0.2
    par.t_cutoff = 800

    par.ωq = 4.47 * 2*pi       # qubit freq
	par.ωc = 7.415 * 2*pi       # cavity freq
	#par.ωd = 7.4143  *2.0*np.pi       # drive freq FOR QB IN GROUND STATE
	#par.ωd = 7.412899  *2.0*np.pi       # drive freq FOR QB IN EXCITED STATE
	par.ωd = 7.4135995 * 2*pi       # drive freq FOR QB IN EXCITED STATE
	par.ah = 0.211 * 2*pi      #qubit anharmonicity
	par.g = 0.0014 * 2*pi     # qubit cavity coupling const
	#par.A = 0.0061 * 2*pi      # drive amplitude
	par.A  = 0.0073 *2*pi      # drive amplitude for about 25 photons in cavity

	if kappa_on
		par.κc = 0.0014 * 2*pi #cavity linewidth
	else
		par.κc = 0.0 * 2*pi #cavity linewidth
	end
	par.κq = 0 *2*pi #qubit linewidth
	par.ϕ = 0
	return par
end

function ket2dm(t, psi)
    return dm(normalize(psi))
end

function ket2dm_qub(t, psi)
    return ptrace(dm(normalize(psi)),2)
end

function mc_evol()
    t_tmp, ρ_avg = timeevolution.mcwf_dynamic(tlist, ψ0, H; fout=ket2dm, maxiters = 1e7)
    for i = 1:Ntrajectories-1
        t_tmp, ρ_temp = timeevolution.mcwf_dynamic(tlist, ψ0, H; fout=ket2dm, maxiters = 1e7)
        ρ_avg = ρ_avg + ρ_temp
        println("$i th iteration complete")
    end

    ρ_avg = ρ_avg./Ntrajectories
    return ρ_avg
end

function mc_evol_qub()
    t_tmp, ρ_avg = timeevolution.mcwf_dynamic(tlist, ψ0, H; fout=ket2dm_qub, maxiters = 1e7)
    for i = 1:Ntrajectories-1
        t_tmp, ρ_temp = timeevolution.mcwf_dynamic(tlist, ψ0, H; fout=ket2dm_qub, maxiters = 1e7)
        ρ_avg = ρ_avg + ρ_temp
        println("$i th iteration complete")
    end

    ρ_avg = ρ_avg./Ntrajectories
    return ρ_avg
end
