using Kinbiont
using SymbolicRegression
using Plots
using Tables
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


using CSV

path_to_data = string("/Fig_3/data/Exp_1/channel_1.csv")
path_to_annotation = string("/Fig_3/data/Exp_1/annotation.csv")
path_to_calib = string("/Fig_3/data//cal_curve_avg.csv")
path_to_results = string("/KinBiont_utilities-main/res/")


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


fit_file =  segmentation_ODE_file(
    "seg_exp_1_calibrated", #label of the experiment
    path_to_data, # path to the folder to analyze
    list_of_models, # ode model to use 
    list_guess, #  param
    1;
    path_to_annotation=path_to_annotation,# path to the annotation of the wells
    detect_number_cpd=false,
    fixed_cpd=false,
    calibration_OD_curve=path_to_calib,
    multiple_scattering_correction=true, # if true uses the given calibration curve to fix the data
    type_of_curve="deriv",
    pt_smooth_derivative= 10,
    type_of_smoothing = "lowess",
    verbose = true,
    write_res =true,
    path_to_results = path_to_results,
    win_size=12, # numebr of the point to generate intial condition
    smoothing =true,
    lb_param_array=list_lb, # lower bound param
    ub_param_array=list_ub, # upper bound param
    maxiters = 200000
)
# changing units of time of gr from 1/min to 1/h
fit_file[2][7,2:end] = fit_file[2][7,2:end] .* 60
fit_file[2][6,2:end] = fit_file[2][6,2:end] .* 60


fit_param =fit_file[2]
index_first_segment = findall(fit_param[end,2:end].==1) .+ 1
index_second_segment  = findall(fit_param[end,2:end].==2) .+ 1
res_first_seg = hcat(fit_param[:,1],fit_param[:,index_first_segment])
res_second_seg = hcat(fit_param[:,1],fit_param[:,index_second_segment])


annotation_test = CSV.File(path_to_annotation,header =false)
names_of_annotation = propertynames(annotation_test)
feature_matrix = hcat(annotation_test[:Column1], annotation_test[:Column3])

# selecting strain S5
index_strain = findall(annotation_test[:Column2].== "S5")
index_cc = findall(annotation_test[:Column3] .> 0.1)
to_keep = intersect(index_strain,index_cc)
feature_matrix = feature_matrix[to_keep,:]
wells_to_use =  feature_matrix[:,1]
index_res =Any
for i in wells_to_use
    if i == wells_to_use[1]
      index_res = findfirst(res_first_seg[:,1:end] .== i)

    else

      index_res = hcat(index_res, findfirst(res_first_seg[:,1:end]  .== i))

    end


end

iii = [ index_res[1,i][2] for i in eachindex( index_res[1,:] )] 

res_first_seg_ML = res_first_seg[:,iii]
res_first_seg_ML = hcat(res_first_seg[:,1],res_first_seg_ML)
# add 0.0 0.0 data for not growing wells
feature_matrix =vcat(feature_matrix,["zero" 0.0])
res_first_seg_ML=hcat(res_first_seg_ML , reduce(vcat,["zero" ,"zero", "zero", 0.0 ,  0.0 ,0.0 ,0.0 ,0.0 ,0.0 ]))



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
    
gr_sy_reg = downstream_symbolic_regression(res_first_seg_ML,
    feature_matrix,
   7;
 #   options = options,
)
 
 
# Example usage
expr_array =string.( gr_sy_reg[1])

x_min = convert(Float64,feature_matrix[1,2])
x_max = convert(Float64,feature_matrix[end,2])
N = 100

results = evaluate_expressions(expr_array, x_min, x_max, N)

scatter( feature_matrix[:,2],res_first_seg_ML[7,2:end],xlabel= "Aminoacid Concentration μM",ylabel = "Growth rate [1/h]",label=[ "Data" nothing])

plot!(results[1],results[2][1],label=[ "Eq. 1" nothing], line=(3,:green,:dash,))
plot!(results[1],results[2][2],label=[ "Eq. 2" nothing], line=(3,:red,))
plot!(results[1] ,results[2][3],label=[ "Eq. 3" nothing],line=(3,:blue,:dashdot,))
plot!(results[1] ,results[2][4],label=[ "Eq. 4" nothing],line=(2,:black,),legend=:bottomright)
 


# selecting strain S6

annotation_test = CSV.File(path_to_annotation,header =false)
names_of_annotation = propertynames(annotation_test)
feature_matrix = hcat(annotation_test[:Column1], annotation_test[:Column3])

index_strain = findall(annotation_test[:Column2].== "S6")
index_cc = findall(annotation_test[:Column3] .> 0.1)
to_keep = intersect(index_strain,index_cc)
feature_matrix = feature_matrix[to_keep,:]
wells_to_use =  feature_matrix[:,1]
index_res =Any
for i in wells_to_use
    if i == wells_to_use[1]
      index_res = findfirst(res_first_seg[:,1:end] .== i)

    else

      index_res = hcat(index_res, findfirst(res_first_seg[:,1:end]  .== i))

    end


end

iii = [ index_res[1,i][2] for i in eachindex( index_res[1,:] )] 

res_first_seg_ML = res_first_seg[:,iii]
res_first_seg_ML = hcat(res_first_seg[:,1],res_first_seg_ML)
# add 0.0 0.0 data
feature_matrix =vcat(feature_matrix,["zero" 0.0])
res_first_seg_ML=hcat(res_first_seg_ML , reduce(vcat,["zero" ,"zero", "zero", 0.0 ,  0.0 ,0.0 ,0.0 ,0.0 ,0.0 ]))



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
    
gr_sy_reg = downstream_symbolic_regression(res_first_seg_ML,
    feature_matrix,
   7;
 #   options = options,
)
 

# Example usage
expr_array =string.( gr_sy_reg[1])

x_min = convert(Float64,feature_matrix[1,2])
x_max = convert(Float64,feature_matrix[end,2])
N = 100

results = evaluate_expressions(expr_array, x_min, x_max, N)

scatter( feature_matrix[:,2],res_first_seg_ML[7,2:end],xlabel= "Aminoacid Concentration [μM]",ylabel = "Growth rate [1/h]",label=[ "Data" nothing])

plot!(results[1],results[2][1],label=[ "Eq. 1" nothing], line=(3,:green,:dash,))
plot!(results[1],results[2][2],label=[ "Eq. 2" nothing], line=(3,:red,))
plot!(results[1] ,results[2][3],label=[ "Eq. 3" nothing],line=(3,:blue,:dashdot,))
plot!(results[1] ,results[2][4],label=[ "Eq. 4" nothing],line=(2,:black,),legend=:bottomright)
 



# N max regression


# selecting strain S5

annotation_test = CSV.File(path_to_annotation,header =false)
names_of_annotation = propertynames(annotation_test)
feature_matrix = hcat(annotation_test[:Column1], annotation_test[:Column3])


index_strain = findall(annotation_test[:Column2].== "S5")
index_cc = findall(annotation_test[:Column3] .> 0.1)
to_keep = intersect(index_strain,index_cc)
feature_matrix = feature_matrix[to_keep,:]
wells_to_use =  feature_matrix[:,1]
index_res =Any
for i in wells_to_use
    if i == wells_to_use[1]
      index_res = findfirst(res_second_seg[:,1:end] .== i)

    else

      index_res = hcat(index_res, findfirst(res_second_seg[:,1:end]  .== i))

    end


end

iii = [ index_res[1,i][2] for i in eachindex( index_res[1,:] )] 

res_second_seg_ML = res_second_seg[:,iii]
res_second_seg_ML = hcat(res_second_seg[:,1],res_second_seg_ML)
# add 0.0 0.0 data
feature_matrix =vcat(feature_matrix,["zero" 0.0])
res_second_seg_ML=hcat(res_second_seg_ML , reduce(vcat,["zero" ,"zero", "zero", 0.0 ,  0.0 ,0.0 ,0.0 ,0.0 ,0.0 ]))


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
    
gr_sy_reg = downstream_symbolic_regression(res_second_seg_ML,
    feature_matrix,
   5;
 #   options = options,
)
# Example usage
expr_array =string.( gr_sy_reg[1])

x_min = convert(Float64,feature_matrix[1,2])
x_max = convert(Float64,feature_matrix[end,2])
N = 100

results = evaluate_expressions(expr_array, x_min, x_max, N)
 

scatter( feature_matrix[:,2],res_second_seg_ML[5,2:end],xlabel= "Aminoacid Concentration [μM]",ylabel ="Total Growth [OD]",label=[ "Data" nothing])

plot!(results[1],results[2][1],label=[ "Eq. 1" nothing], line=(3,:green,:dash,))
plot!(results[1],results[2][2],label=[ "Eq. 2" nothing], line=(3,:red,))
plot!(results[1] ,results[2][3],label=[ "Eq. 3" nothing],line=(3,:blue,:dashdot,))
plot!(results[1] ,results[2][4],label=[ "Eq. 4" nothing],line=(2,:black,),legend=:bottomright)
 


# selecting strain S6

annotation_test = CSV.File(path_to_annotation,header =false)
names_of_annotation = propertynames(annotation_test)
feature_matrix = hcat(annotation_test[:Column1], annotation_test[:Column3])


index_strain = findall(annotation_test[:Column2].== "S6")
index_cc = findall(annotation_test[:Column3] .> 0.1)
to_keep = intersect(index_strain,index_cc)
feature_matrix = feature_matrix[to_keep,:]
wells_to_use =  feature_matrix[:,1]
index_res =Any
for i in wells_to_use
    if i == wells_to_use[1]
      index_res = findfirst(res_second_seg[:,1:end] .== i)

    else

      index_res = hcat(index_res, findfirst(res_second_seg[:,1:end]  .== i))

    end


end

iii = [ index_res[1,i][2] for i in eachindex( index_res[1,:] )] 

res_second_seg_ML = res_second_seg[:,iii]
res_second_seg_ML = hcat(res_second_seg[:,1],res_second_seg_ML)
# add 0.0 0.0 data
feature_matrix =vcat(feature_matrix,["zero" 0.0])
res_second_seg_ML=hcat(res_second_seg_ML , reduce(vcat,["zero" ,"zero", "zero", 0.0 ,  0.0 ,0.0 ,0.0 ,0.0 ,0.0 ]))


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
    
gr_sy_reg = downstream_symbolic_regression(res_second_seg_ML,
    feature_matrix,
   5;
 #   options = options,
)
 

results = evaluate_expressions(expr_array, x_min, x_max, N)
 

scatter( feature_matrix[:,2],res_second_seg_ML[5,2:end],xlabel= "Aminoacid Concentration [μM]",ylabel ="Total Growth [OD]",label=[ "Data" nothing])

plot!(results[1],results[2][1],label=[ "Eq. 1" nothing], line=(3,:green,:dash,))
plot!(results[1],results[2][2],label=[ "Eq. 2" nothing], line=(3,:red,))
plot!(results[1] ,results[2][3],label=[ "Eq. 3" nothing],line=(3,:blue,:dashdot,))
plot!(results[1] ,results[2][4],label=[ "Eq. 4" nothing],line=(2,:black,),legend=:bottomright)
 



