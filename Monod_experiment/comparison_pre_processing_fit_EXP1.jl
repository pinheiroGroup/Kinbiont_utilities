
using Kinbiont


path_to_data = string("/data/Exp_1/channel_1.csv")
path_to_annotation = string("/data/Exp_1/annotation.csv")


path_to_calib = string("/data/cal_curve_avg.csv")
path_to_results = string("/results/")


model = "HPM"

lb_param = [0.00001, 0.000001, 0.01]
ub_param =[0.5,       1.5,     5.5]
param_guess =[0.01, 0.01, 1.02]
    
fit_od = fit_file_ODE(
        "exp_1_no_corrections", #label of the experiment
        path_to_data, # path to the folder to analyze
        model, # string of the used model
        param_guess;
        path_to_annotation=missing,# path to the annotation of the wells
        verbose = true,
        pt_avg = 3, 
        write_res= true,
        path_to_results=path_to_results,
        lb = lb_param,
        ub =ub_param,
        abstol = 0.00000001, 
        maxiters = 20000000,
)

    


    
fit_od = fit_file_ODE(
        "exp_1_BS", #label of the experiment
        path_to_data, # path to the folder to analyze
        model, # string of the used model
        param_guess;
        path_to_annotation=path_to_annotation,# path to the annotation of the wells
        verbose = true,
        pt_avg = 3, 
        write_res= true,
        path_to_results=path_to_results,
        lb = lb_param,
        ub =ub_param,
        abstol = 0.00000001, 
        maxiters = 20000000,
)

    

    
fit_od = fit_file_ODE(
        "exp_1_BS_MS", #label of the experiment
        path_to_data, # path to the folder to analyze
        model, # string of the used model
        param_guess;
        path_to_annotation=path_to_annotation,# path to the annotation of the wells
        verbose = true,
        pt_avg = 3, 
        write_res= true,
        multiple_scattering_correction=true, # if true uses the given calibration curve to fix the data
        calibration_OD_curve=path_to_calib,  #  the path to calibration curve to fix the data
        path_to_results=path_to_results,
        lb = lb_param,
        ub =ub_param,
        abstol = 0.00000001, 
        maxiters = 20000000,
)

    