using Kinbiont
using CSV
using Statistics
using Distributions
using Plots




model_1 = "aHPM"
param_of_ode_1= [0.06, # growth rate
0.01, 
1.01, 
1.5]

# initial contion of the ODE, note sciML requires vectors
n_start =[0.1,0.0]
# starting time of the ODE
tstart =0.0
#ending time of the ODE
tmax = 600.0
# delta t for numerical integration
delta_t=6.0

nrep = 100
# preparing parameters of the fitting
ub_ahpm =[ 0.08 , 0.09 , 1.5  , 5  ]
# lowert bounds of the parameters of the ODE

lb_ahpm =[ 0.01 , 0.01, 0.4 ,0.5 ]

p_guess =[ 0.04 , 0.045 , 0.8  , 1  ]


model_1 = "aHPM"
gr_gt = zeros(nrep)
gr_ODE =zeros(nrep)
gr_data =zeros(nrep)
gr_data_pt2=zeros(nrep)
gr_data_pt3=zeros(nrep)
gr_data_pt4=zeros(nrep)
gr_data_pt5=zeros(nrep)
gr_data_pt6=zeros(nrep)
gr_data_pt7=zeros(nrep)
gr_data_pt8=zeros(nrep)
gr_data_pt9=zeros(nrep)
gr_data_pt10=zeros(nrep)

noise_level = 0.05

noise_levels = [0.0001,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09]
res_final = ["noise","gt","ODE","pt_0","pt_2","pt_3","pt_4","pt_5","pt_6","pt_7","pt_8","pt_9","pt_10"]

path_to_plot = "my_plot_path"

for nn in noise_levels
    noise_level = nn    
    println(string("noise level ", nn))

    for i in 1:nrep 
      println(string("curve ", i))

      param_of_ode_1 =[rand(Uniform(lb_ahpm[i],ub_ahpm[i]),1)[1 ] for i in eachindex(ub_ahpm)]


      sim_1 = ODE_sim(model_1, #string of the model
        n_start, # starting condition
        tstart, # start time of the sim
         tmax, # final time of the sim
            delta_t, # delta t for poisson approx
        param_of_ode_1 # parameters of the ODE model
        
        )
     
     time_sol = reduce(hcat, sim_1.t)
     sol_fin = reduce(hcat, sim_1.u)
     sol_fin = sum(sol_fin,dims=1)
     data_OD =Matrix(vcat(time_sol,sol_fin))
     display(plot(  data_OD[1,:], log.( data_OD[2,:]))  )

      # evaluation gr without noise_unifo
     gr_gt[i] = maximum( specific_gr_evaluation(  Matrix(data_OD), 0))

        # adding noise 

   #  noise_unifom = rand(Normal(0.0,noise_level),length(times_sim))
    noise_unifom = rand(Normal(0,noise_level/2),length(time_sol))

        data_OD[2,:] = data_OD[2,:] .+ noise_unifom
        index_zero =  findall(<(0),  data_OD[2,:] )
    if length(index_zero)>0
        data_OD[2,index_zero] .= 0.1
    end    



     gr_data[i] = maximum( specific_gr_evaluation(  Matrix(data_OD), 0))
     gr_data_pt2[i] = maximum( specific_gr_evaluation(  Matrix(data_OD), 2))
     gr_data_pt3[i] = maximum( specific_gr_evaluation(  Matrix(data_OD), 3))
     gr_data_pt4[i] = maximum( specific_gr_evaluation(  Matrix(data_OD), 4))
     gr_data_pt5[i] = maximum( specific_gr_evaluation(  Matrix(data_OD), 5))
     gr_data_pt6[i] = maximum( specific_gr_evaluation(  Matrix(data_OD), 6))
     gr_data_pt7[i] = maximum( specific_gr_evaluation(  Matrix(data_OD), 7))
     gr_data_pt8[i] = maximum( specific_gr_evaluation(  Matrix(data_OD), 8))
     gr_data_pt9[i] = maximum( specific_gr_evaluation(  Matrix(data_OD), 9))
     gr_data_pt10[i] = maximum( specific_gr_evaluation(  Matrix(data_OD), 10))

        # fitting
    @time fit =  fitting_one_well_ODE_constrained(data_OD, # dataset first row times second row OD
     string(i), # name of the well
     "test_max_gr_detection", #label of the experiment
     "aHPM", # ode model to use
     p_guess; # intial guess 
     lb =lb_ahpm, # lower bound param
     ub = ub_ahpm, # upper bound param
      )
     temp_res = [noise_level,gr_gt[i],gr_ODE[i],gr_data[i],gr_data_pt2[i],gr_data_pt3[i],gr_data_pt4[i],gr_data_pt5[i],gr_data_pt6[i],gr_data_pt7[i],gr_data_pt8[i],gr_data_pt9[i],gr_data_pt10[i]]
     res_final = hcat(res_final,temp_res)
    end
end




