using PyPlot, QuantumOptics, JLD2, Parameters

coupling_type = 3
save_data = true
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

function instate(nq,qub_amp,nc,cav_amp,st_type)
    q_basis = FockBasis(nq)
    c_basis = FockBasis(nc)
    return q_basis, c_basis, ψ0
end

file = jldopen("/home/vamsi/Github/cavity_qubit_system/data_ats 2,2,9.6,0.03330088212805181.jld2", "r")
rho = file["rho"]
par = file["par"]

q_basis = FockBasis(par.nq)
c_basis = FockBasis(par.nc)
include("ham_def.jl")
tlist = file["tlist"]
close(file)

#Cavity Number
exp_nc = real(expect(Npc,rho))
pygui(true)
figure()
plot(exp_nc)
grid("on")
xlabel(L"\mathrm{Time}")
ylabel(L"\mathrm{Photon number}")
title("nq = $(par.nq), nc = $(par.nc)")

#Cavity Population
cav_pop = zeros(length(tlist),par.nc+1)
for i in 0:par.nc
    cav_pop[:,i+1] = real(expect(identityoperator(q_basis) ⊗ dm(fockstate(c_basis,i)),rho))
end
figure()
plot(tlist,cav_pop[:,end-2:end])
xlabel(L"\mathrm{Time}")
ylabel(L"\mathrm{Cavity \, population}")
title("nq = $(par.nq), nc = $(par.nc)")
grid("on")
legend(labels=["$(par.nc-2)" ,"$(par.nc-1)" ,"$(par.nc)"])

#Qubit Population
qub_pop = zeros(length(tlist),par.nq+1)
for i in 0:par.nq
    qub_pop[:,i+1] = real(expect(dm(fockstate(q_basis,i))⊗identityoperator(c_basis),rho))
end
figure()
plot(tlist,qub_pop)
xlabel(L"\mathrm{Time}")
ylabel(L"\mathrm{Qubit \, population}")
title("nq = $(par.nq), nc = $(par.nc)")
grid("on")
legend(labels=["0" ,"1" ,"2", "3" ,"4" ,"5"])

if(save_data == true)
	outfile = "qub_pop_$(par.nc)_photons.txt"
	f = open(outfile, "w")
	println(f,qub_pop)
	close(f)

	outfile = "cav_pop_$(par.nc)_photons.txt"
	f = open(outfile, "w")
	println(f,cav_pop)
	close(f)

	outfile = "parlist_$(par.nc)_photons.txt"
	f = open(outfile, "w")
	println(f,par)
	close(f)

	outfile = "pcavavg_$(par.nc)_photons.txt"
	f = open(outfile, "w")
	println(f,exp_nc)
	close(f)
end
