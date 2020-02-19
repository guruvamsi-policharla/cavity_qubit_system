using Distributed
if (nprocs()==1)
    addprocs(4)
end

@everywhere using QuantumOptics, DifferentialEquations, JLD2, SharedArrays, Parameters
@everywhere include("aux.jl")
@everywhere drive_on = true	 #-- switch for starting in a coherent state, otherwise default is vacuum
@everywhere kappa_on = true  #-- switch decay on or off
@everywhere coupling_type = 3  #-- 1-(Xq)*(Xc), 2-(Npq)*(Npc), 3-(Xq*Xq)*(Xc*Xc)

if (coupling_type == 1)
    println("Transverse chosen")
    @everywhere par = define_params_transverse()
elseif (coupling_type == 2 || coupling_type == 3 || coupling_type == 4 || coupling_type == 5)
    println("Longitudinal chosen")
    @everywhere par = define_params_longitudinal()
    print(par)
else
    error("coupling_type not found")
end

if (drive_on == false)
    println("Initial state is Fock-Coherent")
    par.cav_amp = sqrt(13) #α
    @everywhere q_basis, c_basis, ψ0 = instate(par.nq,par.qub_amp,par.nc,par.cav_amp,"FC")
else
    println("Initial state is Fock-Fock")
    @everywhere q_basis, c_basis, ψ0 = instate(par.nq,par.qub_amp,par.nc,par.cav_amp,"FF")
end

@everywhere include("ham_def.jl")
@everywhere Ntrajectories = 10

ρ_avg = @time @sync @distributed (+) for i = 1:nworkers()
    mc_evol()
end

ρ_avg = ρ_avg./nworkers()
@save "data "*string(par.nq)*","*string(par.nc)*","*string(maximum(tlist))*".jld2" ρ_avg tlist

rmprocs(workers())
println("Successfully removed workers")

#--------------------SAVING FILES-----------------------------------------------
fn = "data "*string(par.nq)*","*string(par.nc)*","*string(maximum(tlist))*".jld2"
jldopen(fn, true, true, true, IOStream; compress=true) do file
    file["rho"] = ρ
    file["tlist"] = tlist
    file["par"] = par
end
close(file)
