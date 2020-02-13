using Distributed
if(nprocs()==1)
    addprocs(2)
end

@everywhere using QuantumOptics, DifferentialEquations, JLD2, SharedArrays, Parameters
@everywhere include("aux.jl")
@everywhere drive_on = True	 #-- switch for starting in a coherent state, otherwise default is vacuum
@everywhere kappa_on = True  #-- switch decay on or off
@everywhere coupling_type = 3  #-- 1-(Xq)*(Xc), 2-(Npq)*(Npc), 3-(Xq*Xq)*(Xc*Xc)

if(coupling_type == 1)
    @everywhere par = define_params_transverse()
elseif(coupling_type == 2 | coupling_type == 3 | coupling_type == 4 | coupling_type == 5)
    @everywhere par = define_params_longitudinal()
end

if(drive_on == flase)
    par.cav_amp = sqrt(13) #α
    @everywhere ψ0 = instate(par.nq,par.qub_amp,par.nc,par.cav_amp,st_type="FC")
else
    @everywhere ψ0 = instate(par.nq,par.qub_amp,par.nc,par.cav_amp,st_type ="FF")
end

include("ham_def.jl")
@everywhere Ntrajectories = 10
if drive_on == true & kappa_on = true
    println("-- DRIVE ON")
	@time @sync @distributed for i = 1:Ntrajectories
	    t_tmp, ψ = timeevolution.mcwf_dynamic(T, ψ0, H; seed=i)
	    exp_nc_average .+= real(expect(Npc, ψ))
	    println(i,"th seed complete")
	end
