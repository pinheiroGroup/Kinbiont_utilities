
using Kinbiont
using CSV,  DataFrames
using Statistics
using Peaks
using ChangePointDetection
using Distributions
using Plots

model_1 = "HPM"
param_of_ode_1= [0.06, # growth rate
0.01, 
1.01, 
]


model_2 = "logistic"
param_of_ode_2= [0.08,
2.0, 
]
# initial contion of the ODE, note sciML requires vectors
n_start =[0.1,0.0]
# starting time of the ODE
tstart =0.0
#ending time of the ODE
tmax = 600.0
# delta t for numerical integration
delta_t=6.0

nrep =100

times_gt = zeros(nrep)
times_detected =zeros(nrep)


noise_levels = [0.0001,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09]
res_final = ["noise","gt","detected","peak_prominence","#"]



times_gt = zeros(nrep)
times_detected =zeros(nrep)
times_detected_p =zeros(nrep)

for i in 1:nrep 

    t_cp_s = rand(Uniform(150,450),1)[1]
    times_gt[i]=t_cp_s
    sim_1 = ODE_sim(model_1, #string of the model
     n_start, # starting condition
     tstart, # start time of the sim
     t_cp_s, # final time of the sim
     delta_t, # delta t for poisson approx
     param_of_ode_1 # parameters of the ODE model
     )
     
     time_sol = reduce(hcat, sim_1.t)
     sol_fin = reduce(hcat, sim_1.u)
     sol_fin = sum(sol_fin,dims=1)
     # t chage point
     sim_2 = ODE_sim(model_2, #string of the model
     [sol_fin[end]], # starting condition
     t_cp_s, # start time of the sim
     tmax, # final time of the sim
     delta_t, # delta t for poisson approx
     param_of_ode_2 # parameters of the ODE model
     )
 
     sol_2 =reduce(hcat,sim_2.u)
 
     times_sim =vcat(sim_1.t,sim_2.t)
     sol_sim =hcat(sol_fin,sol_2)
     data_OD =Matrix(vcat(transpose(times_sim),sol_sim))
     data_OD_noise = copy(data_OD)
     for nn in noise_levels
        noise_level = nn    
        println(string("noise level ", nn))
    

        noise_unifom = rand(Normal(0,noise_level/2),length(data_OD[2,:]))
        data_OD_noise[2,:] = data_OD[2,:] .+ noise_unifom
 

        cdp =cpd_local_detection(data_OD_noise,
             1; type_of_curve ="original", method="thr_scan",number_of_bin=300)

        cdp_p =cpd_local_detection(data_OD_noise,
             1; type_of_curve ="original",size_win =2, method="peaks_prominence")
   


        res_temp_1 = [noise_level,t_cp_s,cdp[2][1],cdp_p[2][1],1]
        res_final =hcat(res_final,res_temp_1) 


     end
 end 

scatter(res_final[2,2:end],res_final[3,2:end])
scatter!(res_final[2,2:end],res_final[4,2:end])



