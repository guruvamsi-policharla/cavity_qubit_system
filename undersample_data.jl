using QuantumOptics, JLD2

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
fn = "data 5,9,20.0.jld2"
file = jldopen(fn, "r")
rho = file["rho"]
tlist = file["tlist"]
println(length(rho))
println(length(tlist))
rho = rho[1:4:end]
tlist = tlist[1:4:end]
println(length(rho))
println(length(tlist))



jldopen("undersampled_"*fn, true, true, true, IOStream; compress=true) do file
    file["rho"] = rho
    file["tlist"] = tlist
    file["par"] = par
end
close(file)
