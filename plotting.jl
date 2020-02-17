using QuantumOptics
using JLD2
using PyPlot
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

file = jldopen("/home/vamsi/Github/cavity_qubit_system/data 5,25,1000.0.jld2", "r")
ρ = file["ρ"]

#Cavity Number
exp_nc_1 = real(expect(Npc,ρ))
pygui(true)
figure()
plot(tlist,exp_nc_1)
grid("on")
xlabel(L"\mathrm{Time}")
ylabel(L"\mathrm{Photon number}")

#Qubit Population
qub_pop = zeros(length(tlist),par.nq+1)
for i in 0:par.nq
    qub_pop[:,i+1] = qub_0 = real(expect(dm(fockstate(q_basis,i))⊗identityoperator(c_basis),ρ))
end
figure()
plot(tlist,qub_pop)
