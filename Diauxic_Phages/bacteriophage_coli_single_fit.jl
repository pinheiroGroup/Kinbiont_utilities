using Kinbiont
using CSV
using Tables
using StatsBase
using Plots

################################
list_files  = ["10^8cfu"]

blank_value= mean([0.0965,	0.0943,	0.0974, 0.0963,	0.0947,	0.0972,0.0964,	0.0943,	0.0957])


ub_dhpm_d =[ 0.1 , 0.7, 0.1  ,0.1]
lb_dhpm_d =[ 0.000000001 , 0.0000000010, 0.000000001  ,0.000000000001]
p_guess  = lb_dhpm_d .+ (ub_dhpm_d.- lb_dhpm_d)./2

model = "HPM_3_death"

# Please change the paths 
path_to_data = string("E:/Lavoro/Kinbiont_utilities-main/Kinbiont_utilities-main/Diauxic_Phages///data/coli_phages/10^8cfu.csv")
path_to_results = string("//")


# number of segment

data_matrix  = CSV.read(path_to_data, Tables.matrix)
data_OD = permutedims( [data_matrix[:,1] data_matrix[:,12]] )
data_OD[2,:] = data_OD[2,:] .- blank_value
scatter(data_OD[1,:],data_OD[2,:], label="Data", xlabel="Time [h]", ylabel="OD [arb. units]", legend=:topleft)




results_ODE_fit = Kinbiont.fitting_one_well_ODE_constrained(
    data_OD, 
    "test",
    "Coli_Phages",
    model,
    p_guess;
   lb = lb_dhpm_d,
   ub = ub_dhpm_d
)



scatter(data_OD[1,:],data_OD[2,:], label="Data", xlabel="Time [h]", ylabel="OD [arb. units]", legend=:topleft,color =:black)
plot!(results_ODE_fit[4],results_ODE_fit[3],label=[ "Fit" nothing],line=4,xlabel="Time [h]", ylabel="OD [arb. units]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,alpha = 0.5,color =:blue)
savefig(string("fitting_coli_phages_examples.svg"))