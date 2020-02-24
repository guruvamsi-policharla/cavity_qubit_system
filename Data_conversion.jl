using QuantumOptics, JLD2

file = jldopen("/home/vamsi/Github/cavity_qubit_system/Data/Transverse/24-2/13 Photons/data 5,26,1000.0,0.03330088212805181.jld2", "r")
rho = file["rho"]
tlist = file["tlist"]

outfile = "data 5,26,1000.0,0.03330088212805181.txt"
f = open(outfile, "w")

for i in 1:length(rho)
    println(f,rho[i].data)
    print(f,"\n")
end
close(f)

outfile = "tlist_small.txt"
f = open(outfile, "w")
for i in 1:length(tlist)
    println(f,tlist)
end
close(f)
