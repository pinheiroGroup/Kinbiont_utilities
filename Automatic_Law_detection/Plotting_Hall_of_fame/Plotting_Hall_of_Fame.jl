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


#plot!(results[1],results[2][1],label=[ "Eq. I" nothing], line=(3,:green,:dash,),xlabel= "Amino acid Concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11)
plot(results[1],results[2][2],label=[ "Eq. II" nothing], line=(6,:blue,:dashdot,),xlabel= "Amino acid concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,alpha = 0.5)
plot!(results[1] ,results[2][3],label=[ "Eq. III" nothing],line=(6,:black),xlabel= "Amino acid concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11)
plot!(results[1],results[2][4],label=[ "Eq. IV" nothing], line=(6,:blue,:dash,),xlabel= "Amino acid concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,alpha = 0.5)
scatter!( concentration_array,parse.(Float64,first_segment_well[7,:]),xlabel= "Amino acid concentration [μM]",ylabel = "Growth rate [1/h]",label=[ "Data" nothing],size = (600,500),markersize = 6,color= :green,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,legend=:bottomright)
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


#plot!(results[1],results[2][1],label=[ "Eq. I" nothing], line=(3,:green,:dash,),xlabel= "Amino acid Concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11)
plot(results[1],results[2][2],label=[ "Eq. X" nothing], line=(6,:blue,:dashdot,),xlabel= "Amino acid concentration [μM]",ylabel = "Total growth [OD]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,alpha = 0.5)
plot!(results[1] ,results[2][3],label=[ "Eq. XI" nothing],line=(6,:black),xlabel= "Amino acid concentration [μM]",ylabel = "Total Growth [OD]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11)
plot!(results[1],results[2][4],label=[ "Eq. XII" nothing], line=(6,:blue,:dash,),xlabel= "Amino acid concentration [μM]",ylabel = "Total growth [OD]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,alpha = 0.5)
scatter!( concentration_array,parse.(Float64,second_segment_well[5,:]),xlabel= "Amino acid concentration [μM]",ylabel = "Total growth [OD]",label=[ "Data" nothing],size = (600,500),markersize = 6,color= :red,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,legend=:bottomright)
savefig(string("../Monod_experiment/Nmax_S5.svg"))

 # Plotting results for CHL dose response


 
# loading results and annotation with the concentration of antibiotic
results_of_fit = CSV.read("../Dose_response_sym_regression/Hall_of_fame/results_full.csv", Tables.matrix)
annotation = CSV.read("../Dose_response_sym_regression/Data/annotation_full.csv", Tables.matrix)


# finding wells with CHL

index_not_used_wells  = findall(annotation[:,2].== "X" .|| annotation[:,2].== "b")
non_blank_wells = setdiff(1:1:length(annotation[:,1]),index_not_used_wells)
annotation_test_no_blank= annotation[non_blank_wells,:]
index_rif_wells  = annotation_test_no_blank[findall( annotation_test_no_blank[:,5].== "Rifampicin" ),1]
index_df_wells  =  annotation_test_no_blank[findall( annotation_test_no_blank[:,5].== "DF" ),1]
index_chl_wells  =  annotation_test_no_blank[findall( annotation_test_no_blank[:,5].== "Chloramphenicol" ),1]


index_to_use_df = [ findall(  index_df_wells[i,1] .==annotation_test_no_blank[:,1] ) for i in 1:length(index_df_wells[:,1])] 
index_to_use_df = reduce(vcat,index_to_use_df)

index_to_use_rif =  [ findall( index_rif_wells[i,1] .== annotation_test_no_blank[:,1] ) for i in 1:length(index_rif_wells[:,1])] 
index_to_use_rif =  reduce(vcat,index_to_use_rif)
index_to_use_rif =  vcat(index_to_use_df,index_to_use_rif)

index_to_use_chl =   [ findall( index_chl_wells[i,1] .== annotation_test_no_blank[:,1] ) for i in 1:length(index_chl_wells[:,1])] 
index_to_use_chl =  reduce(vcat,index_to_use_chl)
index_to_use_chl =  vcat(index_to_use_df,index_to_use_chl)



feature_matrix_chl =hcat(annotation_test_no_blank[index_to_use_chl,1], annotation_test_no_blank[index_to_use_chl,6])
feature_matrix_rif =hcat(annotation_test_no_blank[index_to_use_rif,1], annotation_test_no_blank[index_to_use_rif,6])


# selecting the growth rate matrix for rif 
# selectin index of results matrix to use for rif and chl

index_to_use_chl_results = [ findall(  feature_matrix_chl[i,1] .== results_of_fit[2,:]) for i in 1:length( feature_matrix_chl[:,1])] 
index_to_use_chl_results =  reduce(vcat,index_to_use_chl_results)

index_to_use_rif_results = [ findall(  feature_matrix_rif[i,1] .== results_of_fit[2,:]) for i in 1:length( feature_matrix_chl[:,1])] 
index_to_use_rif_results =  reduce(vcat,index_to_use_rif_results)



results_matrix_chl =hcat(results_of_fit[:,1] ,results_of_fit[:,index_to_use_chl_results])
results_matrix_rif =hcat(results_of_fit[:,1] ,results_of_fit[:,index_to_use_rif_results])


concetrations_chl = feature_matrix_chl[:,2]
concetrations_rif = feature_matrix_rif[:,2]



# plotting chl hall of fame
max_conc = convert(Float64,maximum(concetrations_chl))
min_conc = 0.00


Hall_of_Fame = CSV.read("../Dose_response_sym_regression/Hall_of_fame/hall_of_fame_chl.csv", Tables.matrix)
expr_array =string.( Hall_of_Fame[:,3])

# number of point to plot the function
N = 100

results = evaluate_expressions(expr_array, min_conc, max_conc, N)


#plot!(results[1],results[2][1],label=[ "Eq. I" nothing], line=(3,:green,:dash,),xlabel= "Amino acid Concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11)
plot(results[1],results[2][2],label=[ "Eq. XX" nothing], line=(6,:blue,:dashdot,),xlabel= "Chloramphenicol concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,alpha = 0.5)
plot!(results[1] ,results[2][3],label=[ "Eq. XXI" nothing],line=(6,:black),xlabel= "Chloramphenicol concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11)
plot!(results[1],results[2][4],label=[ "Eq. XXII" nothing], line=(6,:blue,:dash,),xlabel= "Chloramphenicol concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,alpha = 0.5)
scatter!( concetrations_chl,parse.(Float64,results_matrix_chl[6,2:end]),xlabel= "Chloramphenicol concentration [μM]",ylabel = "Growth rate [1/h]",label=[ "Data" nothing],size = (600,500),markersize = 6,color= :green,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,legend=:topright)
savefig(string("../Dose_response_sym_regression/GR_VS_CHL.svg"))




# plotting rif hall of fame
max_conc = convert(Float64,maximum(concetrations_rif))
min_conc = 0.00


Hall_of_Fame = CSV.read("../Dose_response_sym_regression/Hall_of_fame/hall_of_fame_rif.csv", Tables.matrix)
expr_array =string.( Hall_of_Fame[:,3])

# number of point to plot the function
N = 100

results = evaluate_expressions(expr_array, min_conc, max_conc, N)


#plot!(results[1],results[2][1],label=[ "Eq. I" nothing], line=(3,:green,:dash,),xlabel= "Amino acid Concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11)
plot(results[1],results[2][2],label=[ "Eq. XX" nothing], line=(6,:blue,:dashdot,),xlabel= "Rifampicin concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,alpha = 0.5)
plot!(results[1] ,results[2][3],label=[ "Eq. XXI" nothing],line=(6,:black),xlabel= "Rifampicin concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11)
plot!(results[1],results[2][4],label=[ "Eq. XXII" nothing], line=(6,:blue,:dash,),xlabel= "Rifampicin concentration [μM]",ylabel = "Growth rate [1/h]",size = (600,500),markersize = 6,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,alpha = 0.5)
scatter!( concetrations_rif,parse.(Float64,results_matrix_rif[6,2:end]),xlabel= "Rifampicin concentration [μM]",ylabel = "Growth rate [1/h]",label=[ "Data" nothing],size = (600,500),markersize = 6,color= :green,tickfontsize = 20,labelfontsize = 20,legendfontsize =11,legend=:topright)
savefig(string("../Dose_response_sym_regression/GR_VS_RIF.svg"))
