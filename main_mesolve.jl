using QuantumOptics, DifferentialEquations, JLD2, SharedArrays, Parameters
include("aux.jl")
drive_on = true	 #-- switch for starting in a coherent state, otherwise default is vacuum
kappa_on = true  #-- switch decay on or off
coupling_type = 3  #-- 1-(Xq)*(Xc), 2-(Npq)*(Npc), 3-(Xq*Xq)*(Xc*Xc)

if (coupling_type == 1)
    println("Transverse chosen")
    par = define_params_transverse()
elseif (coupling_type == 2 || coupling_type == 3 || coupling_type == 4 || coupling_type == 5)
    println("Longitudinal chosen")
    par = define_params_longitudinal()
    print(par)
else
    error("coupling_type not found")
end

if (drive_on == false)
    println("Initial state is Fock-Coherent")
    par.cav_amp = sqrt(13) #α
    q_basis, c_basis, ψ0 = instate(par.nq,par.qub_amp,par.nc,par.cav_amp,"FC")
else
    println("Initial state is Fock-Fock")
    q_basis, c_basis, ψ0 = instate(par.nq,par.qub_amp,par.nc,par.cav_amp,"FF")
end

include("ham_def.jl")

@time tout, ρ = timeevolution.master_dynamic(tlist, ψ0, H; maxiters = 1e7)

#--------------------SAVING FILES-----------------------------------------------
fn = "data "*string(par.nq)*","*string(par.nc)*","*string(maximum(tlist))*".jld2"
jldopen(fn, true, true, true, IOStream; compress=true) do file
    file["rho"] = ρ
    file["tlist"] = tlist
    file["par"] = par
end
close(file)
