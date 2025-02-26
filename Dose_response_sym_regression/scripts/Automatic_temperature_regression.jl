# loading data form  "Modeling of Bacterial Growth as a Function of Temperature"
using Kinbiont
using Plots
using Statistics
using Distributions
using Tables
using SymbolicRegression
using CSV


Data_matrix = CSV.File("/Users/fabrizio.angaroni/Documents/KinBiont_utilities-main/temperature_dose_response/plot-data.csv")
names_of_cols = propertynames(Data_matrix)
res_of_fitting = Data_matrix[names_of_cols[2]]
growth_rate =  vcat("gr",Data_matrix[names_of_cols[2]])
labels_names = [ string(i) for i in 1:length( Data_matrix[names_of_cols[2]])]
labels_names =  vcat("labels",labels_names)
results_for_ML = hcat(labels_names,growth_rate)
results_for_ML = permutedims(hcat(labels_names,results_for_ML))

# add 0.0 0.0 data for not growing wells
feature_matrix = Data_matrix[names_of_cols[1]]
feature_matrix =hcat(labels_names[2:end], feature_matrix)



# define SymbolicRegression options
options = SymbolicRegression.Options(;
 binary_operators=[+,/,-,* ]  ,
 constraints = [(+) =>(-1,1),(-) =>(-1,1),(/) =>(-1,1),(*) =>(-1,1)],
 unary_operators=[exp,square],
 #constraints=nothing,
 elementwise_loss=nothing,
 loss_function=nothing,
 tournament_selection_n=14, #1 sampled from every tournament_selection_n per mutation
 tournament_selection_p=0.86,
 topn=14, #samples to return per population
 complexity_of_operators=nothing,
 complexity_of_constants=nothing,
 complexity_of_variables=nothing,
 parsimony=0.0000000005,
 dimensional_constraint_penalty=nothing,
 nested_constraints = [exp => [exp => 0],square => [square => 0],exp => [square => 0], square => [exp => 0]],

 alpha=0.00000001,
 optimizer_probability = 1.0,
 maxsize=12,
 optimizer_iterations = 100,
 should_optimize_constants=true,
 maxdepth=nothing,
 seed =50,
)



# regression on growth rate
# using as variable 1/t to try to obtaion 
#feature_matrix[:,2] = exp.(.- 1 ./feature_matrix[:,2])
gr_sy_reg = downstream_symbolic_regression(results_for_ML,
    feature_matrix,
   3;
   options = options,
)
 
gr_sy_reg[1]

scatter(feature_matrix[:,2],results_for_ML[3,2:end])
hline!(unique(gr_sy_reg[3][:,1]),label=[ "Eq. 1" nothing], line=(3,:green,:dash,),xlabel ="exp(-1/Temperature)",ylabel="Growth rate [1/h]")
plot!(unique(convert.(Float64,feature_matrix[gr_sy_reg[4],2])) ,unique(gr_sy_reg[3][:,2]),label=[ "Eq. 2" nothing], line=(3,:red,))
plot!(unique(convert.(Float64,feature_matrix[gr_sy_reg[4],2])) ,unique(gr_sy_reg[3][:,3]),label=[ "Eq. 3" nothing],line=(3,:blue,:dashdot,))
plot!(unique(convert.(Float64,feature_matrix[gr_sy_reg[4],2])) ,unique(gr_sy_reg[3][:,4]),label=[ "Eq. 4" nothing],line=(2,:black,))
plot!(unique(convert.(Float64,feature_matrix[gr_sy_reg[4],2])) ,unique(gr_sy_reg[3][:,5]),label=[ "Eq. 5" nothing],line=(2,:green,))
plot!(unique(convert.(Float64,feature_matrix[gr_sy_reg[4],2])) ,unique(gr_sy_reg[3][:,6]),label=[ "Eq. 6" nothing],line=(2,:violet,))
plot!(unique(convert.(Float64,feature_matrix[gr_sy_reg[4],2])) ,unique(gr_sy_reg[3][:,7]),label=[ "Eq. 7" nothing],line=(2,:gray,))
plot!(unique(convert.(Float64,feature_matrix[gr_sy_reg[4],2])) ,unique(gr_sy_reg[3][:,8]),label=[ "Eq. 8" nothing],line=(2,:black,))
plot!(unique(convert.(Float64,feature_matrix[gr_sy_reg[4],2])) ,unique(gr_sy_reg[3][:,9]),label=[ "Eq. 9" nothing],line=(2,:black,))

gr_sy_reg[1]

