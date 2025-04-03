using Kinbiont
using Statistics
using Distributions
using Tables
using SymbolicRegression
using CSV
using Random
using Plots

#################### 

# setting the path to the data and the annotation file
path_to_calibration = "../Data/cal_curve_avg.csv"
path_to_data =        "../Data/data_H.csv"
path_to_annotation =  "../Data/annotation_short.csv"
path_to_plot = "../Fit_plots/"
# creating the path to save the plots

mkpath(path_to_plot)
# setting parameters of the models

model1 = "HPM_exp"

lb_param1 = [0.00001, 0.000001]
ub_param1 =[3.5,       10.5]
param_guess1 =[0.1, 0.01]
    
model2 = "logistic"

lb_param2 = [0.00001, 0.000001]
ub_param2 =[4.5,       5.5    ]
param_guess2 =[0.1, 1.01]
    
# creating the array of models, param_guess, lb_param, ub_param for Kinbiont

list_of_models = [model1,model2]
list_guess=  [param_guess1,param_guess2]
list_lb=  [lb_param1,lb_param2]
list_ub=  [ub_param1,ub_param2]

# calling the segmentation function, we fix the number of cp to 1 to separate the growth rate from the stationary phase

Kinbiont_seg_analysis =  segmentation_ODE_file(
  "fit_antibio_response", #label of the experiment
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
  pt_smooth_derivative= 12,
  verbose = true,
  write_res =false,
  type_of_smoothing="rolling_avg",
  win_size=12, # 
  smoothing =true,
  lb_param_array=list_lb, # lower bound param
  ub_param_array=list_ub, # upper bound param
  maxiters = 200000
)
## save the fits as svg
# a kinbiont fit result data structure 
# contains the fits in the third element
Kinbiont_fits =Kinbiont_seg_analysis[3];
# the data in the fourth element
Kinbiont_data =Kinbiont_seg_analysis[4];
# the results matrix of the fits in the second element
Kinbiont_results_matrix =Kinbiont_seg_analysis[2];
# the change point list in the fifth element
Kinbiont_cp_intervals =Kinbiont_seg_analysis[5];
# the label of the wells is in the second row of the results matrix
well_names = unique(Kinbiont_seg_analysis[2][2,2:end]);


# we plot each well with a for loop
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
         ylabel="arb. units",
         label=["Data " nothing],
         color=:black,
         guidefontsize=15,
         tickfontsize=15,
         markersize=8,
         legendfontsize=15,
         size=(600,500),
     ),
     )

     display(
      Plots.plot!(
        x_fit_temp,
          y_fit_temp,
          xlabel="Time [h]",
          ylabel="OD [arb. Units]",
          label=["fit " nothing],
          line = 8,
          alpha = 0.8,
          color=:blue,
          guidefontsize=15,
          tickfontsize=15,
          legendfontsize=15,
          size=(600,500),
      ),
      )


     display( Plots.vline!(
         [temp_cp],
         c=:black,
         label=[string("Change point") nothing],
         guidefontsize=15,
         linesize = 8,
         tickfontsize=15,
         legendfontsize=15,
         legend=:bottomright,
         size=(600,500),
         ) )
         savefig(string(path_to_plot,  "segmented_fit_", well_name, ".svg"))
end
## 

index_first_segment = findall(Kinbiont_seg_analysis[2][9,:].==1)

results_matrix =hcat(Kinbiont_seg_analysis[2][:,1] ,Kinbiont_seg_analysis[2][:,index_first_segment])


# write results matrix
CSV.write("../Data/results_full.csv",Tables.table(results_matrix))


# we inport the annotation containing the concentrations of the drugs and the type of drug
# selecting the wells drug free and chl  in the annotation file that contains also the concentrations


path_to_annotation_full=  "../Data/annotation_full.csv"
annotation_test = CSV.read(path_to_annotation_full,header =false,Tables.matrix)


index_not_used_wells  = findall(annotation_test[:,2].== "X" .|| annotation_test[:,2].== "b")
non_blank_wells = setdiff(1:1:length(annotation_test[:,1]),index_not_used_wells)
annotation_test_no_blank= annotation_test[non_blank_wells,:]
index_df_wells  =  annotation_test_no_blank[findall( annotation_test_no_blank[:,5].== "DF" ),1]
index_chl_wells  =  annotation_test_no_blank[findall( annotation_test_no_blank[:,5].== "Chloramphenicol" ),1]


index_to_use_df = [ findall(  index_df_wells[i,1] .==annotation_test_no_blank[:,1] ) for i in 1:length(index_df_wells[:,1])] 
index_to_use_df = reduce(vcat,index_to_use_df)


index_to_use_chl =   [ findall( index_chl_wells[i,1] .== annotation_test_no_blank[:,1] ) for i in 1:length(index_chl_wells[:,1])] 
index_to_use_chl =  reduce(vcat,index_to_use_chl)
index_to_use_chl =  vcat(index_to_use_df,index_to_use_chl)



feature_matrix_chl =hcat(annotation_test_no_blank[index_to_use_chl,1], annotation_test_no_blank[index_to_use_chl,6])


# selectin index of results matrix to use for  chl

index_to_use_chl_results = [ findall(  feature_matrix_chl[i,1] .== results_matrix[2,:]) for i in 1:length( feature_matrix_chl[:,1])] 
index_to_use_chl_results =  reduce(vcat,index_to_use_chl_results)



results_matrix_chl =hcat(results_matrix[:,1] ,results_matrix[:,index_to_use_chl_results])





scatter(feature_matrix_chl[:,2],results_matrix_chl[6,2:end],xlabel= "CHL Concentration [μM]",ylabel = "Gr [1/h]",label=[ "Data" nothing],size = (600,500),markersize = 6,color= :red,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,legend=:topright)

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


# regression on growth rate chl


gr_chl = downstream_symbolic_regression(results_matrix_chl,
feature_matrix_chl,
   6;
 #   options = options,
);
 
gr_chl[1]