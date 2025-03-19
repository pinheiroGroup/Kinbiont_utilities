using Kinbiont
using Plots
using Tables
using CSV
using Tables

path_to_data = string("../data/Exp_1/channel_1.csv")
path_to_annotation = string("../data/Exp_1/annotation.csv")
path_to_calib = string("../data/cal_curve_avg.csv")

data_matrix  = CSV.read(path_to_data, Tables.matrix)
data_OD = permutedims( [data_matrix[:,1] data_matrix[:,12]] )
data_OD[2,:] = data_OD[2,:] .-0.08
scatter(data_OD[1,:],data_OD[2,:], label="Data", xlabel="Time [h]", ylabel="OD [arb. units]", legend=:topleft)

model1 = "HPM_exp"

lb_param1 = [0.00001, 0.000001]
ub_param1 =[0.5,       1.5]
param_guess1 =[0.01, 0.01]
    
model2 = "logistic"

lb_param2 = [0.00001, 0.000001]
ub_param2 =[0.5,       2.5,    ]
param_guess2 =[0.01, 1.01]
    

list_of_models = [model1,model2]
list_guess=  [param_guess1,param_guess2]
list_lb=  [lb_param1,lb_param2]
list_ub=  [ub_param1,ub_param2]



seg_fitting = Kinbiont.segmentation_ODE(
   data_OD, # dataset first row times second row OD
    "test", # name of the well
    "", #label of the experiment
    list_of_models, #  models to use
    list_guess,
    1;
    calibration_OD_curve=path_to_calib,
    multiple_scattering_correction=true, # if true uses the given calibration curve to fix the data
    type_of_curve="deriv",
    pt_smooth_derivative= 10,
    type_of_smoothing = "lowess",
    verbose = true,
    win_size=12, # numebr of the point to generate intial condition
    smoothing =true,
    lb_param_array=list_lb, # lower bound param
    ub_param_array=list_ub, # upper bound param
    maxiters = 2000000,
)



scatter(data_OD[1,:],data_OD[2,:], label="Data", xlabel="Time [h]", ylabel="OD [arb. units]", legend=:topleft,color =:black)
vline!(seg_fitting[5],label=[ "Change point" nothing],line=3,color=:black)
plot!(seg_fitting[4],seg_fitting[3],label=[ "Fit" nothing],line=4,xlabel="Time [h]", ylabel="OD [arb. units]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,alpha = 0.5,color =:blue)
savefig(string("../fitting.svg"))