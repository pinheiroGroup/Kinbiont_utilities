using Kinbiont
using CSV
using Tables
using Plots
################################

ub_exp = [0.3]
lb_exp = [00.00]

ub_logistic = [0.1, 0.9]
lb_logistic = [0.0001, 0.000001]

ub_alogistic = [0.1, 0.9, 10.0]
lb_alogistic = [0.0001, 0.000001, 0.001]

ub_hpm = [0.3, 0.8, 0.9]
lb_hpm = [0.0001, 0.000001, 0.001]


ub_hpm_exp = [0.1, 1.0]
lb_hpm_exp = [0.0001, 0.0000001]



list_of_models = ["HPM", "HPM_exp", "logistic", "exponential"]
list_ub_param = [ub_hpm, ub_hpm_exp, ub_logistic, ub_exp]
list_lb_param = [lb_hpm, lb_hpm_exp, lb_logistic, lb_exp]
list_guess = list_lb_param .+ (list_ub_param .- list_lb_param)./2

data_1 = CSV.File("E:/Lavoro/Kinbiont_utilities-main/Kinbiont_utilities-main/Diauxic_Phages/data/diauxic/data_modified_1.csv", transpose=true, missingstring=nothing)





names_of_cols = propertynames(data_1)
data_temp = data_1[names_of_cols[1]]
times_data = convert(Vector{Float64}, data_temp)



data_OD = Matrix(transpose(hcat(times_data, data_1[names_of_cols[40]])))
scatter(data_OD[1,:],data_OD[2,:], label="Data", xlabel="Time [h]", ylabel="OD [arb. units]", legend=:topleft)

# a well with 2 segment
seg_fitting = segmentation_ODE(
            data_OD, # dataset first row times second row OD
            string(r), # name of the well
            "data_1", #label of the experiment
            list_of_models, # ode models to use
            list_guess,
            2;
            lb_param_array=list_lb_param, # lower bound param
            ub_param_array=list_ub_param, # upper bound param
            detect_number_cpd=false,
            fixed_cpd=false,
            type_of_detection="slinding_win",
            type_of_curve="deriv",
            pt_avg=3, # number of the point to generate intial condition
            smoothing=false, # the smoothing is done or not?            win_size=10, # numebr of the point to generate intial condition
            pt_smooth_derivative=6,
            win_size = 10
)




scatter(data_OD[1,:],data_OD[2,:], label="Data", xlabel="Time [h]", ylabel="OD [arb. units]", legend=:topleft,color =:black)
vline!(seg_fitting[5],label=[ "Change point" nothing],line=3,color=:black)
plot!(seg_fitting[4],seg_fitting[3],label=[ "Fit" nothing],line=4,xlabel="Time [h]", ylabel="OD [arb. units]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,alpha = 0.5,color =:blue)
savefig(string("../fitting_diauxic_examples_2_seg_big.svg"))