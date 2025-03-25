using Kinbiont
using Statistics
using Distributions
using Tables
using SymbolicRegression
using CSV
using Random
using Plots

#################### 

path_to_calibration = "../Data/cal_curve_avg.csv"
path_to_data =        "../Data/data_channel_1.csv"
path_to_annotation =  "../Data/Annotation_CHL_dosage.csv"
path_to_plot = "../Fit_plots/Annotation_CHL_dosage.csv"
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


Kinbiont_seg_analysis =  segmentation_ODE_file(
  "fit_chl_response", #label of the experiment
  path_to_data, # path to the folder to analyze
  list_of_models, # ode model to use 
  list_guess, #  param
  1;
  path_to_annotation=path_to_annotation,# path to the annotation of the wells
  detect_number_cpd=false,
  fixed_cpd=false,
  calibration_OD_curve=path_to_calibration,
  multiple_scattering_correction=true, # if true uses the given calibration curve to fix the data
  type_of_curve="deriv",
  pt_smooth_derivative= 10,
  verbose = true,
  write_res =false,
  win_size=12, # numebr of the point to generate intial condition
  smoothing =true,
  lb_param_array=list_lb, # lower bound param
  ub_param_array=list_ub, # upper bound param
  maxiters = 200000
)
## save the fits as svg
Kinbiont_fits =Kinbiont_seg_analysis[3]

Kinbiont_results_matrix =Kinbiont_seg_analysis[2]
Kinbiont_data =Kinbiont_seg_analysis[4]
Kinbiont_cp_intervals =Kinbiont_seg_analysis[5]
well_names = unique(Kinbiont_seg_analysis[2][2,2:end])



for i in eachindex(Kinbiont_fits)
   fit_temp =Kinbiont_fits[i]
   data_temp =Kinbiont_data[i]
   temp_cp =Kinbiont_cp_intervals[i]
   y_fit_temp =fit_temp[:,2]
   x_fit_temp =fit_temp[:,1]
   y_data_temp =data_temp[2,:]
   x_data_temp =data_temp[1,:]
   well_name = well_names[i]

   display(
     Plots.scatter(
         x_data_temp,
         y_data_temp,
         xlabel="Time",
         ylabel="Arb. Units",
         label=["Data " nothing],
         markersize=4,
         color=:black,
         guidefontsize=15,
         tickfontsize=15,
         legendfontsize=15,
         size=(400,300),
     ),
     )

     display(
      Plots.plot(
         x_data_temp,
          y_fit_temp,
          xlabel="Time",
          ylabel="Arb. Units",
          label=["Data " nothing],
          markersize=4,
          linesize = 4,
          color=:red,
          guidefontsize=15,
          tickfontsize=15,
          legendfontsize=15,
          size=(400,300),
      ),
      )


     display( Plots.vline!(
         [temp_cp],
         c=:black,
         label=[string("Change point") nothing],
         guidefontsize=guidefontsize,
         tickfontsize=15,
         legendfontsize=15,
         size=(400,300),
         ) )
         savefig(string(path_to_plot,  "segmented_fit_", well_name, ".svg"))
end
## 

index_first_segment = findall(Kinbiont_seg_analysis[2][9,:].==1)

results_matrix =hcat(Kinbiont_seg_analysis[2][:,1] ,Kinbiont_seg_analysis[2][:,index_first_segment])


# write results matrix
CSV.write("../Data/results_matrix_CHL_dosage.csv",Tables.table(results_matrix))

# selecting chl


path_to_annotation_full= path_to_annotation
annotation_test = CSV.File(path_to_annotation_full,header =false)

names_of_annotation = propertynames(annotation_test)
feature_matrix = hcat(annotation_test[:Column1], annotation_test[:Column5])

index_not_used_wells  = findall(annotation_test[:Column2].== "X" .|| annotation_test[:Column2].== "b")
index_to_use = setdiff(1:1:length(annotation_test[:Column2]),index_not_used_wells)



res_of_fitting = results_matrix
# add 0.0 0.0 data for not growing wells
feature_matrix =feature_matrix[index_to_use,:]



# define SymbolicRegression options
options = SymbolicRegression.Options(;
 binary_operators=[+,  /, * , -]  ,
 unary_operators=[],
 constraints=nothing,
 elementwise_loss=nothing,
 loss_function=nothing,
 tournament_selection_n=12, #1 sampled from every tournament_selection_n per mutation
 tournament_selection_p=0.86,
 topn=12, #samples to return per population
 complexity_of_operators=nothing,
 complexity_of_constants=nothing,
 complexity_of_variables=nothing,
 parsimony=0.05,
 dimensional_constraint_penalty=nothing,
 alpha=0.100000,
 maxsize=10,
 seed =1234,
 maxdepth=nothing,
 

)

# regression on growth rate


gr_sy_reg = downstream_symbolic_regression(res_of_fitting,
    feature_matrix,
   6;
 #   options = options,
)
 
gr_sy_reg[1]
