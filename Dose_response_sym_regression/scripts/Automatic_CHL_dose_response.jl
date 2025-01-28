using Kinbiont
using KinbiontPlots
using Plots
using Statistics
using Distributions
using Tables
using SymbolicRegression
using CSV
# adding function to plots please refer to correct file



#################### 

path_to_plots =  "/Users/fabrizio.angaroni/Documents/KinBiont_utilities-main/Dose_response_sym_regression/Plots/"
path_to_res = "/Users/fabrizio.angaroni/Documents/KinBiont_utilities-main/Dose_response_sym_regression/Results/"
path_to_calibration = "/Users/fabrizio.angaroni/Documents/KinBiont_utilities-main/Dose_response_sym_regression/Data/cal_curve_avg.csv"
path_to_data = "/Users/fabrizio.angaroni/Documents/KinBiont_utilities-main/Dose_response_sym_regression/Data/data_channel_1.csv"
path_to_annotation =  "/Users/fabrizio.angaroni/Documents/KinBiont_utilities-main/Dose_response_sym_regression/Data/annotation_channel_1_media_M9+glucose.csv"


    


Kinbiont_seg_analysis =  segment_gr_analysis_file(
    path_to_data,
    string("CHL_RIF_dose_response"); 
    path_to_annotation = path_to_annotation,
    path_to_results=path_to_res,
    write_res=true,
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



plot_fit_of_file(
    Kinbiont_seg_analysis;
        path_to_plot=path_to_plots, # path where to save Plots
        display_plots=true,# display plots in julia or not
        save_plots=true, # save the plot or not
        x_size=500,
        pt_smoothing_derivative = 10,
        y_size =750,
        guidefontsize=20,
        tickfontsize=18,
        legendfontsize=12,
)
index_first_segment = findall(Kinbiont_seg_analysis[2][11,:].==1)

results_matrix = Kinbiont_seg_analysis[2][:,index_first_segment]
results_matrix = hcat( Kinbiont_seg_analysis[2][:,1],results_matrix)





# selecting chl

path_to_annotation_full= "/Users/fabrizio.angaroni/Documents/KinBiont_utilities-main/Dose_response_sym_regression/Data/Annotation_CHL_dosage.csv"
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
 maxdepth=nothing,
 

)

# regression on growth rate
    
gr_sy_reg = downstream_symbolic_regression(res_of_fitting,
    feature_matrix,
   3;
 #   options = options,
)
 
gr_sy_reg[1]



using SymPy

function evaluate_expressions(expr_array::Vector{String}, x_min::Float64, x_max::Float64, N::Int)
    # Define the symbolic variable
    x1 = symbols("x1")
    
    # Convert character array into symbolic expressions
    expressions = [sympy.sympify(expr) for expr in expr_array]
    
    # Generate N points between x_min and x_max
    x_vals = range(x_min, x_max, length=N)
    
    # Evaluate each expression at each point in x_vals
    results = [Float64[subs(expr, x1 => x) for x in x_vals] for expr in expressions]
    
    return x_vals,results
end

# Example usage
expr_array =string.( gr_sy_reg[1])

x_min = convert(Float64,feature_matrix[1,2])
x_max = convert(Float64,feature_matrix[end,2])
N = 100

results = evaluate_expressions(expr_array, x_min, x_max, N)

scatter( feature_matrix[:,2],res_of_fitting[3,2:end],xlabel= "Amino Acid concentration μM",ylabel = "Growth rate [1/Min]",label=[ "Data" nothing])

plot!(results[1],results[2][1],label=[ "Eq. 1" nothing], line=(3,:green,:dash,))
plot!(results[1],results[2][2],label=[ "Eq. 2" nothing], line=(3,:red,))
plot!(results[1] ,results[2][3],label=[ "Eq. 3" nothing],line=(3,:blue,:dashdot,))
plot!(results[1] ,results[2][4],label=[ "Eq. 4" nothing],line=(2,:black,))
plot!(results[1] ,results[2][5],label=[ "Eq. 5" nothing],line=(2,:black,))

using SymbolicUtils

using SymbolicUtils
using Symbolics

function find_matching_expressions(expr_array::Vector{String}, pattern_str::String)
    # Define the symbolic variable
    @variables x1

    # Parse the pattern string into a symbolic expression
    pattern_expr = Meta.parse(pattern_str)
    pattern = Symbolics.build_function(pattern_expr)

    # Initialize an array to hold indices of matching expressions
    matching_indices = []

    # Iterate over the expression strings
    for (i, expr_str) in enumerate(expr_array)
        # Parse the expression string into a symbolic expression
        expr = Meta.parse(expr_str)
        symbolic_expr = Symbolics.build_function(expr)
        SymbolicUtils.simplify(expr)
        # Check if the symbolic expression matches the pattern
        if SymbolicUtils.ismatch(pattern, symbolic_expr)
            push!(matching_indices, i)
        end
    end

    return matching_indices
end

using Symbolics: match
desired_form = 1 / (1 + x1)

match_result = match(desired_form, SymbolicUtils.simplify(expr))


SymbolicUtils.ismatch(pattern,  SymbolicUtils.simplify(expr))

pattern_str = "1 / (1 + x1)"
matching_indices = find_matching_expressions(expr_array, pattern_str)
println("Expressions matching the pattern found at indices: ", matching_indices)
symbolic_expr = SymbolicUtils.substitute(expr, Dict(:x1 => x1))
