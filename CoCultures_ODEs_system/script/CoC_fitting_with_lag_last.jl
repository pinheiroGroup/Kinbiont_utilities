using Kinbiont
using CSV
using DataFrames
using Plots
using Statistics
# data from "Construction and Modeling of a Coculture Microplate for Real-Time Measurement of Microbial Interactions"
# strain annotation 
A0 = [1,11,21,31,41,51,13,33,53,7,17,37,47,57,38,48,58]
LP = [12,32,52,5,15,25,35,45,55,27,8,9,19,29,10,20,30]
LB = [4,14,24,34,44,54,16,36,56,18,28,39,49,59,40,50,60]

# load data
data = CSV.read("E:/Lavoro/Kinbiont_utilities-main/CoCultures_ODEs_system/data/DATA_Gut_Microbiome_CoC.csv", DataFrame)
# convet time column from in elapsed time frome the start in hours, the delta time in the data is 15 minutes

delta_time = 15/60
time_seq = 0:delta_time:(size(data,1)-1)*delta_time
time_seq = [ time_seq[i] for i in 1: length(time_seq)]
data[!,:Time] = time_seq 

index_tot = 2:1:size(data,2)
A0_index = A0.+1
LP_index = LP.+1
LB_index = LB.+1
blanS1_index = setdiff(index_tot, [A0_index; LP_index; LB_index])

### evaluating blank values by doing the mean of the blank columns
data_blanS1_values = Matrix(data[:,blanS1_index])
blanS1_value = mean(reduce(vcat,eachrow(data_blanS1_values)))
blanS1_value = 0.0

##########################################################
## Selecting two wells with interactions 
##########################################################




# A0 and LP
A0_inct =[13,33,53,17].+1
LP_inct =[13,33,53,17].+2

data_A0_inct = data[:,A0_inct] .- blanS1_value
data_LP_inct = data[:,LP_inct] .-  blanS1_value

# mean of the replicates
data_A0_inct_mean = mean(Matrix(data_A0_inct),dims=2)

data_LP_inct_mean = mean(Matrix(data_LP_inct),dims=2)
scatter(data_A0_inct_mean)
scatter!(data_LP_inct_mean)

# Formatting data for Kinbiont
data_to_fit = time_seq
data_to_fit = permutedims( [time_seq data_A0_inct_mean data_LP_inct_mean])



# CoC model

function CoC_dynamics(du, u, param, t)
    B_1, B_2,G = u
   mu_1,mu_2, interaction_1_2, interaction_2_1,lag_1,lag_2,K_A0_G,K_LP_G= param


  
    # flux of growth with limiting

    VA0_GL = G/(K_A0_G + G)
    VLP_GL = G/(K_LP_G + G)
    lag_factor1 = (t/(t+lag_1))
    lag_factor2 = (t/(t+lag_2))

    
    
    du[1] =lag_factor1*( mu_1 *( VA0_GL)* B_1  + interaction_1_2 *  du[2]  )# dB_1 # biomass A0
    du[2] =lag_factor2* (mu_2 * ( VLP_GL) * B_2   + interaction_2_1 * du[1] )# dB_2 # biomass t LP
    du[3] = -lag_factor1* VA0_GL * B_1 - lag_factor2*VLP_GL * B_2  # dG_L_dt glucose

end

u0 = [mean(data_to_fit[2,1:3]),mean(data_to_fit[3,1:3]),2*0.0055]  # Initial conditions NOTE THE DOUBLE OF THE INITIAL GLUCOSE



param = [655.362,63.8116,10.243, -0.19,0.946,24,87,5]  #  K_A0_G, K_LP_G, K_A0_S1, K_LP_S2,leak_1,leak_2, mu_max_1,mu_max_22 
ub_p = [1000.0,1000.0,100.0, 100.0,100,100,100,100]
lb_p = [000.0,00.0,-100.0, -100.0,0.0,0.0,0.0,0.0]
# Initial condition of the states

# parameters 


# fitting

fit = fit_ODEs_System(
    data_to_fit,
    "test",  # Label for dataset
    CoC_dynamics,  # Custom ODE function
    param,  # Initial parameter guess
    u0 ; # Initial conditions
set_equattion_to_fit= [1,2],
ub=ub_p,
lb=lb_p,
maxiters=100000)


fit[2]
sol = fit[3]



Plots.scatter(data_to_fit[1,:],data_to_fit[2,:], xlabel="Time [h]", ylabel="OD [Arb. Units]", label=["Data AO" nothing],color=:black,markersize =2 ,size = (300,300))
Plots.plot!(data_to_fit[1,:],sol[1,:], xlabel="Time [h]", ylabel="OD [Arb. Units]",label=["fit AO" nothing],color=:blue,markersize =4 ,size = (600,500),legendposition = :topleft,linewidth=4,linestyle=:dashdot,tickfontsize = 20,labelfontsize = 20,legendfontsize =11)

Plots.scatter!(data_to_fit[1,:],data_to_fit[3,:], xlabel="Time [h]", ylabel="OD [Arb. Units]", label=["Data LP" nothing],color=:green,markersize =2 ,size = (300,300),marker=:diamond)
Plots.plot!(data_to_fit[1,:],sol[2,:], xlabel="Time [h]", ylabel="OD [Arb. Units]",label=["fit LP" nothing],color=:green,markersize =4 ,size = (600,500),legendposition = :topleft,linewidth=4,linestyle=:dashdot,tickfontsize = 20,labelfontsize = 20,legendfontsize =11)

savefig(string("E:/Lavoro/Kinbiont_utilities-main/CoC_fitting_double.svg"))
