using KinBiont
using CSV


nl_model = ["NL_Richards"]
p_guess = [[0.8,1.0,0.1, 40.0]]
lb_nl =[[0.01,0.0001,0.000001,00.001]]
ub_nl =[[3.01,5.01,3.000001,300.01]]


path_to_annotation = string("/data_isolate_chem_mixture_clean/annotation.csv")
path_to_data0 = "/data_isolate_chem_mixture_clean/clean_data/"
list_of_data_chem = readdir(path_to_data0)

for plate in list_of_data_chem

    println(plate)
    path_to_data = string(path_to_data0,plate)
    temp_name = convert(String,split(plate, ".")[1])


    path_to_results = string("/res/Richards//", temp_name, "//")





    fit_nl = fit_NL_model_selection_file(
        string("Richards_", temp_name), #label of the experiment
    path_to_data    , # path to the folder to analyze
    nl_model, # ode model to use
    p_guess;# initial guess param
    lb_param_array =lb_nl, # lower bound param
    ub_param_array = ub_nl, # upper bound param
    path_to_annotation = path_to_annotation,# path to the annotation of the wells
#    multistart = true,
  #  n_restart = 10,
    maxiters = 200000,
    abstol = 0.000001,
    write_res = true,
    path_to_results = path_to_results
  
)



end