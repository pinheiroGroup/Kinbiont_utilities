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
data = CSV.read("../data/DATA_Gut_Microbiome_CoC.csv", DataFrame)
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
##########################################################
## Selecting two wells without interactions  
##########################################################

function CMR_GLU_dynamics(du, u, param, t)
    G,  B_1 = u
    K_A0_G, mu_max_1 = param


  
    # flux of growth with limiting

    VA0_GL = G/(K_A0_G + G)


    du[1] = - VA0_GL * B_1  # dG_L_dt glucose on the left

    
    du[2] = mu_max_1 *( VA0_GL )* B_1   # dB_L_dt # biomass in the left A0
  
end


# A0 
A0_inct = [1,21,41].+1

data_A0_inct = data[:,A0_inct] .- blanS1_value

# mean of the replicates
data_A0_inct_mean = mean(Matrix(data_A0_inct),dims=2)

scatter(data_A0_inct_mean)


# Formatting data for Kinbiont
data_to_fit = permutedims( [time_seq data_A0_inct_mean])



# formulatig the differential model for Kinbiont
# single species model




u0 = [0.0055, 0.01]  # Initial conditions
param = [0.001, 0.1]  #  
ub_p = [500.5, 1000.0]
lb_p = [0.0, 0.0]



fit_AO = fit_ODEs_System(
    data_to_fit,
    "A0",  # Label for dataset
    CMR_GLU_dynamics,  # Custom ODE function
    param,  # Initial parameter guess
    u0 ; # Initial conditions
    set_of_equations_to_fit= [2],
ub=ub_p,
lb=lb_p
)


fit_AO[2]
sol= fit_AO[3]
scatter(data_to_fit[2,:])
plot!(sol[2,:])

# fitting LP


# A0 
LP_inct = [5,25,45].+1

data_LP_inct = data[:,LP_inct] .- blanS1_value

# mean of the replicates
data_LP_inct_mean = mean(Matrix(data_LP_inct),dims=2)

scatter(data_A0_inct_mean)


# Formatting data for Kinbiont
data_to_fit = permutedims( [time_seq data_LP_inct_mean])







fit_LP = fit_ODEs_System(
    data_to_fit,
    "A0",  # Label for dataset
    CMR_GLU_dynamics,  # Custom ODE function
    param,  # Initial parameter guess
    u0 ; # Initial conditions
    set_of_equations_to_fit= [2],
ub=ub_p,
lb=lb_p
)

fit_LP[2]
sol= fit_LP[3]
scatter(data_to_fit[2,:])

plot!(sol[2,:])
##########################################################
## Selecting two wells with interactions 
##########################################################


# taking the fixed parameters from th RESULTS of the single species model
mu_max_1 = fit_AO[2][2,4]
mu_max_2 = fit_LP[2][2,4]
K_A0_G = fit_AO[2][2,3]
K_LP_G = fit_LP[2][2,3]




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
    interaction_1_2, interaction_2_1 = param


  
    # flux of growth with limiting

    VA0_GL = G/(K_A0_G + G)
    VLP_GL = G/(K_LP_G + G)

    
    
    du[1] = mu_max_1 *( VA0_GL)* B_1  + interaction_1_2 *  du[2]  # dB_1 # biomass A0
    du[2] = mu_max_2 * ( VLP_GL) * B_2   + interaction_2_1 * du[1] # dB_2 # biomass t LP
    du[3] = - VA0_GL * B_1 - VLP_GL * B_2  # dG_L_dt glucose

end

u0 = [2*0.0055, 0.01,0.01]  # Initial conditions NOTE THE DOUBLE OF THE INITIAL GLUCOSE
param = [0.1, 0.1]  #  K_A0_G, K_LP_G, K_A0_S1, K_LP_S2,leak_1,leak_2, mu_max_1,mu_max_22 
ub_p = [100.0, 100.0]
lb_p = [-100.0, -100.0]
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
#maxiters=10000
)

fit[2]
sol = fit[3]
scatter(data_A0_inct_mean)
scatter!(data_LP_inct_mean)
plot!(sol[1,:])
plot!(sol[2,:])