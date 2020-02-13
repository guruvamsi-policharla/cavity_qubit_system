using QuantumOptics
using JLD2
using PyPlot

function basis_gen(Nq,Nc)
    if mod(Nq,2) == 0
        q_basis = SpinBasis((Nq-1) // 2)
    else
        q_basis = SpinBasis((Nq-1)/2 // 1)
    end

    c_basis = FockBasis(Nc)

    return q_basis, c_basis
end

nq = 5
nc = 9
q_basis, c_basis = basis_gen(nq,nc)

aq = sigmam(q_basis) ⊗ identityoperator(c_basis)    # qubit annihilation op
ac = identityoperator(q_basis) ⊗ destroy(c_basis)    # cavity annihilation op

Iqc = identityoperator(q_basis) ⊗ identityoperator(c_basis) # Identity operator

Xq = (aq + dagger(aq)); Yq = -1*im*(aq - dagger(aq)); Zq = (Iqc - 2*dagger(aq)*aq); Npq = dagger(aq)*aq;
Xc = (ac + dagger(ac)); Yc = -1*im*(ac - dagger(ac)); Zc = (Iqc - 2*dagger(ac)*ac); Npc = dagger(ac)*ac;
T = 0:0.01:100

#@load "/home/vamsi/Github/cavity_qubit_system/data 5,9,100.0.jld2" exp_nc_average
file = jldopen("/home/vamsi/Github/cavity_qubit_system/Data/data 5,49,100.0,700.jld2", "r")
exp_nc_1 = file["_exp_nc"]
file = jldopen("/home/vamsi/Github/cavity_qubit_system/Data/data 5,49,100.0,800.jld2", "r")
exp_nc_2 = file["_exp_nc"]
pygui(true)
figure()
plot(T, exp_nc_1, label="700")
plot(T, exp_nc_2, label="800")
grid("on")
xlabel(L"\mathrm{Time}")
ylabel(L"\mathrm{Photon number}")
title("5,49")
legend(loc="upper left")
