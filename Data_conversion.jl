using QuantumOptics, JLD2
for i in 1:4
    folder = "/home/vamsi/Github/cavity_qubit_system/Data/MonteCarlo_qub/Transverse/"
    file = jldopen(folder*"data 5,26,1000.0_$(i).jld2", "r")
    rho = file["rho_qub"]
    tlist = file["tlist"]
    par = file["par"]
    outfile = "data 5,26,1000.0_$(i).txt"
    f = open(outfile, "w")

    for i in 1:length(rho)
        println(f,rho[i].data)
        print(f,"\n")
    end
    close(f)

    outfile = "tlist26.txt"
    f = open(outfile, "w")
    for i in 1:length(tlist)
        println(f,tlist)
    end
    close(f)

    outfile = "par26.txt"
    f = open(outfile, "w")
    println(f,par)
    close(f)
end
