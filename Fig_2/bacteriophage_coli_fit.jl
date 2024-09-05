using Kinbiont
using CSV
using Tables
using StatsBase

################################
list_files  = ["10^8cfu"]

ub_exp =[ 0.6 ]
lb_exp =[ -00.1 ]

ub_logistic =[ 0.1 , 1.501]
lb_logistic =[ 0.0001 , 0.001]


ub_hpm =[ 0.1 , 2.0 , 5.001  ]
lb_hpm =[ 0.0001 , 0.000001, 0.001  ]


ub_hpm_exp = [ 0.1 , 0.7   ]
lb_hpm_exp = [ 0.0001 , 0.0000001  ]


ub_dhpm_d =[ 0.1 , 0.7, 0.1  ,0.1]
lb_dhpm_d =[ 0.000000001 , 0.0000000010, 0.000000001  ,0.000000000001]


list_of_models = ["exponential","HPM","HPM_exp","logistic","HPM_3_death"]
list_lb_param =[ub_exp,ub_hpm,ub_hpm_exp,ub_logistic,ub_dhpm_d]
list_ub_param =[lb_exp,lb_hpm,lb_hpm_exp,lb_logistic,lb_dhpm_d]
list_guess = list_lb_param .+ (list_ub_param .- list_lb_param)./2

# Please change the paths 
path_to_data = string("../KinBiont_utilities-main/Fig_2/data/coli_phages/10^8cfu.csv")
path_to_results = string("//")


# number of segment
segment_number = [0, 1, 2, 3 ]


for nseg in segment_number


    fit_file = segmentation_ODE_file(
    string("segmented_ode_", nseg), #label of the experiment
    path_to_data, # path to the folder to analyze
    list_of_models, # ode model to use
    list_guess,
    nseg;
    lb_param_array=list_ub_param, # lower bound param
    ub_param_array=list_lb_param, # upper bound param
    detect_number_cpd=false,
    fixed_cpd=false,
    type_of_detection="sliding_win",
    type_of_curve="deriv",
    do_blank_subtraction="avg_blank",
    correct_negative="remove",
    smoothing=false, # the smoothing is done or not?
    path_to_results=path_to_results,
    win_size = 8, # numebr of the point to generate intial condition
    pt_smooth_derivative=6,
    write_res= true,
    type_of_smoothing="lowess",
    verbose=true,
    blank_value= mean([0.0965,	0.0943,	0.0974, 0.0963,	0.0947,	0.0972,0.0964,	0.0943,	0.0957]),
    blank_array = [0.0965,	0.0943,	0.0974, 0.0963,	0.0947,	0.0972,0.0964,	0.0943,	0.0957],)
    # writing results
    #writing fits
    CSV.write(string(path_to_results,"/fits_coli_phages_",nseg,".csv"),Tables.table( [fit_file[3][k] for k in eachindex(fit_file[3])]))
    # writing data pre-processed
    CSV.write(string(path_to_results,"/data_coli_phages_",nseg,".csv"),Tables.table( [fit_file[4][k] for k in eachindex(fit_file[4])]))
    # writing AIC
    CSV.write(string(path_to_results,"/AIC_coli_phages_",nseg,".csv"),Tables.table( [fit_file[end][k] for k in eachindex(fit_file[end])]))
    # writing cp
    CSV.write(string(path_to_results,"/CP_coli_phages_",nseg,".csv"),Tables.table( [fit_file[5][k] for k in eachindex(fit_file[5])]))


end


 