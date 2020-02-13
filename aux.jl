@everywhere function instate(nq,qub_amp=0,nc,cav_amp=0,st_type="FF")
	#Nq = levels in the qubit
	#qub_amp = amplitude vector of qubit state
	#Nc = levels in the cavity
	#st_type = type of state to be prepared. Fock-Fock by default
	#ca_amp = amplitude vector of cavity state/a scalar alpha if coherence state
    q_basis = FockBasis(nq)
    c_basis = FockBasis(nc)
    if(st_type=="FF")
        ψ0 = Ket(q_basis, qub_amp) ⊗ fockstate(c_basis,cav_amp)
    elseif(st_type == "FC")
        ψ0 = Ket(q_basis, qub_amp) ⊗ coherentstate(c_basis,cav_amp)
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
            κ::Float64 = 0
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
    par.nq = 4
    par.nc = 40
    par.qub_amp=[0,1,0,0]
    par.t = 1000
    par.ts = 20
    par.t_cutoff = 800

    par.ωq = 4.47 * 2*pi       # qubit freq
	par.ωc = 7.415 * 2*pi       # cavity freq
	par.ωd = 7.4205 * 2*pi       # drive freq FOR QB IN EXCITED STATE
	#par.ωd = 7.4218*2.0*np.pi       # drive freq FOR QB IN GROUND STATE
	par.ah = 0.211 * 2*pi      #qubit anharmonicity
	par.g = 0.163 * 2*pi     # qubit cavity coupling const
	par.A1 = sqrt(2) * (0.1*1.4/27) * 2*pi      # drive amplitude

	if kappa_on:
		par.κc = 0.0014 * 2*pi #cavity linewidth
	else:
		par.κc = 0.0 * 2*pi #cavity linewidth

	par.κq = 0 *2*pi #qubit linewidth
	par.ϕ = 0
	print(par)
	return par
end

function define_params_longitudinal()
    par = param_list()
    par.nq = 4
    par.nc = 40
    par.qub_amp=[0,1,0,0]
    par.t = 100
    par.ts = 20
    par.t_cutoff = 80

    par.ωq = 4.47 * 2*pi       # qubit freq
	par.ωc = 7.415 * 2*pi       # cavity freq
	#par.ωd = 7.4143  *2.0*np.pi       # drive freq FOR QB IN GROUND STATE
	#par.ωd = 7.412899  *2.0*np.pi       # drive freq FOR QB IN EXCITED STATE
	par.ωd = 7.4135995 * 2*pi       # drive freq FOR QB IN EXCITED STATE
	par.ah = 0.211 * 2*pi      #qubit anharmonicity
	par.g = 0.0014 * 2*pi     # qubit cavity coupling const
	par.A1 = 0.0061 * 2*pi      # drive amplitude

	if kappa_on:
		par.κc = 0.0014 * 2*pi #cavity linewidth
	else:
		par.κc = 0.0 * 2*pi #cavity linewidth

	par.κq = 0 *2*pi #qubit linewidth
	par.ϕ = 0
	print(par)
	return par
end
