using QuantumOptics, JLD2

file = jldopen("/home/vamsi/Github/cavity_qubit_system/Data/data 5,9,10.0.jld2", "r")
ρ = file["ρ"]
tlist = file["tlist"]

outfile = "small_data.txt"
f = open(outfile, "w")

for i in 1:length(ρ)
    println(f,ρ[i].data)
    print(f,"\n")
end
close(f)

outfile = "tlist_small.txt"
f = open(outfile, "w")
for i in 1:length(tlist)
    println(f,tlist)
end
close(f)
