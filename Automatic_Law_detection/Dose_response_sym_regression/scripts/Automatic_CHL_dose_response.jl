using Kinbiont
using Statistics
using Distributions
using Tables
using SymbolicRegression
using CSV
using Random


#################### 

path_to_calibration = "../Data/cal_curve_avg.csv"
path_to_data =        "../Data/data_channel_1.csv"
path_to_annotation =  "../Data/Annotation_CHL_dosage.csv"

Kinbiont_seg_analysis =  segment_gr_analysis_file(
    path_to_data,
    string("CHL_dose_response"); 
    path_to_annotation = path_to_annotation,
   #type_of_smoothing="lowess",
    type_of_detection="slinding_win",
    multiple_scattering_correction=true,
    calibration_OD_curve=path_to_calibration,
    type_of_curve="deriv",
    win_size=16,
    pt_avg =14,
   # thr_negative=0.01,
    pt_smoothing_derivative=10,
    n_max_change_points = 1, 
    correct_negative="remove", 
   #thr_lowess=0.01,
)




index_first_segment = findall(Kinbiont_seg_analysis[2][11,:].==1)

results_matrix = Kinbiont_seg_analysis[2][:,index_first_segment]
results_matrix = hcat( Kinbiont_seg_analysis[2][:,1],results_matrix)


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
   3;
 #   options = options,
)
 
gr_sy_reg[1]

