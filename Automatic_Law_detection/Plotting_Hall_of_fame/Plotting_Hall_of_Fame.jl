using CSV
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
    
    return x_vals,result

end

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
 
