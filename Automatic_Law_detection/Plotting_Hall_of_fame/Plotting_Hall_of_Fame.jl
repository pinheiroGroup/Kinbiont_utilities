using CSV
using Plots
using Tables
using DataFrames
using SymPy
using Tables
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
####################################
# Monod experiment data
#####################################

results_of_fit = CSV.read("../Monod_experiment/Hall_of_Fame/fit_results.csv", Tables.matrix)
annotation = CSV.read("../Monod_experiment/Data/Exp_1/annotation.csv", Tables.matrix)
# selecting wells with S5
names_S5 = annotation[findall(x -> occursin("S5", x), annotation[:,2]),1]
# find S5 wells int the results of the fit
S5_wells = results_of_fit[:,findall(x -> x in names_S5, results_of_fit[2,:])]
first_segment_well = S5_wells[:,findall(x -> x == "1", S5_wells[end,:])]
concentration_array = annotation[findall(x -> occursin("S5", x), annotation[:,2]),3]
# find the concentration associated to the wells
max_conc = maximum(annotation[findall(x -> occursin("S5", x), annotation[:,2]),3])
min_conc = 0.001
# Plotting GR S5 
Hall_of_Fame = CSV.read("../Monod_experiment/Hall_of_Fame/hall_of_fame_GR_S5.csv", Tables.matrix)
expr_array =string.( Hall_of_Fame[:,3])

# number of point to plot the function
N = 100

results = evaluate_expressions(expr_array, min_conc, max_conc, N)


#plot!(results[1],results[2][1],label=[ "Eq. I" nothing], line=(3,:green,:dash,),xlabel= "Aminoacid Concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11)
plot(results[1],results[2][2],label=[ "Eq. II" nothing], line=(6,:blue,:dashdot,),xlabel= "Aminoacid Concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,alpha = 0.5)
plot!(results[1] ,results[2][3],label=[ "Eq. III" nothing],line=(6,:black),xlabel= "Aminoacid Concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11)
plot!(results[1],results[2][4],label=[ "Eq. IV" nothing], line=(6,:blue,:dash,),xlabel= "Aminoacid Concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,alpha = 0.5)
scatter!( concentration_array,parse.(Float64,first_segment_well[7,:]),xlabel= "Aminoacid Concentration [μM]",ylabel = "Growth rate [1/h]",label=[ "Data" nothing],size = (600,500),markersize = 6,color= :green,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,legend=:bottomright)
savefig(string("../Monod_experiment/GR_S5.svg"))

 
# Plotting N_max S5

names_S5 = annotation[findall(x -> occursin("S5", x), annotation[:,2]),1]
# find S5 wells int the results of the fit
S5_wells = results_of_fit[:,findall(x -> x in names_S5, results_of_fit[2,:])]
second_segment_well = S5_wells[:,findall(x -> x == "2", S5_wells[end,:])]
concentration_array = annotation[findall(x -> occursin("S5", x), annotation[:,2]),3]
# find the concentration associated to the wells
max_conc = maximum(annotation[findall(x -> occursin("S5", x), annotation[:,2]),3])
min_conc = 0.001
# Plotting GR S5 
Hall_of_Fame = CSV.read("../Monod_experiment/Hall_of_Fame/hall_of_fame_N_max_S5.csv", Tables.matrix)
expr_array =string.( Hall_of_Fame[:,3])

# number of point to plot the function
N = 100

results = evaluate_expressions(expr_array, min_conc, max_conc, N)


#plot!(results[1],results[2][1],label=[ "Eq. I" nothing], line=(3,:green,:dash,),xlabel= "Aminoacid Concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11)
plot(results[1],results[2][2],label=[ "Eq. X" nothing], line=(6,:blue,:dashdot,),xlabel= "Aminoacid Concentration [μM]",ylabel = "Total Growth [OD]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,alpha = 0.5)
plot!(results[1] ,results[2][3],label=[ "Eq. XI" nothing],line=(6,:black),xlabel= "Aminoacid Concentration [μM]",ylabel = "Total Growth [OD]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11)
plot!(results[1],results[2][4],label=[ "Eq. XII" nothing], line=(6,:blue,:dash,),xlabel= "Aminoacid Concentration [μM]",ylabel = "Total Growth [OD]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,alpha = 0.5)
scatter!( concentration_array,parse.(Float64,second_segment_well[5,:]),xlabel= "Aminoacid Concentration [μM]",ylabel = "Total Growth [OD]",label=[ "Data" nothing],size = (600,500),markersize = 6,color= :red,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,legend=:bottomright)
savefig(string("../Monod_experiment/Nmax_S5.svg"))

 # Plotting results for CHL dose response


 

results_of_fit = CSV.read("../Dose_response_sym_regression/Hall_of_fame/results_matrix_CHL_dosage.csv", Tables.matrix)
annotation = CSV.read("../Dose_response_sym_regression/Data/Annotation_CHL_dosage.csv", Tables.matrix)


# findin wells with CHL



index_CHL_wells = setdiff(eachindex(annotation[:,5]),findall(ismissing,annotation[:,5]))
concetrations = annotation[index_CHL_wells,5]
scatter( concetrations,parse.(Float64,results_of_fit[3,2:end]),xlabel= "CHL conc [μM]",ylabel = "Growth rate [1/h]",label=[ "Data" nothing],size = (600,500),markersize = 6,color= :red,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,legend=:bottomright)

max_conc = convert(Float64,maximum(concetrations))
min_conc = 0.00


Hall_of_Fame = CSV.read("../Dose_response_sym_regression/Hall_of_fame/hall_of_fame_CHL_dose_response.csv", Tables.matrix)
expr_array =string.( Hall_of_Fame[:,3])

# number of point to plot the function
N = 100

results = evaluate_expressions(expr_array, min_conc, max_conc, N)


#plot!(results[1],results[2][1],label=[ "Eq. I" nothing], line=(3,:green,:dash,),xlabel= "Aminoacid Concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11)
plot(results[1],results[2][2],label=[ "Eq. XX" nothing], line=(6,:blue,:dashdot,),xlabel= "Chloramphenicol Concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,alpha = 0.5)
plot!(results[1] ,results[2][3],label=[ "Eq. XXI" nothing],line=(6,:black),xlabel= "Chloramphenicol Concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11)
plot!(results[1],results[2][4],label=[ "Eq. XXII" nothing], line=(6,:blue,:dash,),xlabel= "Chloramphenicol Concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,alpha = 0.5)
scatter!( concetrations,parse.(Float64,results_of_fit[3,2:end]),xlabel= "Chloramphenicol Concentration [μM]",ylabel = "Growth rate [1/h]",label=[ "Data" nothing],size = (600,500),markersize = 6,color= :green,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,legend=:topright)
savefig(string("../Dose_response_sym_regression/GR_VS_CHL.svg"))
