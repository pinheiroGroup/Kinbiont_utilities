
# Modeling and kinetic study of bio-ethanol production from soy protein concentrate by-product
using Kinbiont
using Plots
using StatsBase
# Data from 
# Time points
time_data = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 28.0, 32.0, 36.0, 40.0, 44.0]

# Experimental data
cells = exp.([7.51, 7.51, 7.59, 7.66, 7.75, 7.89, 7.96, 8.08, 8.13, 8.14, 8.14, 8.15, 8.15, 8.15, 8.16, 8.16, 8.16, 8.16])
exp_TS = [170.20, 165.62, 161.15, 159.55, 150.67, 140.66, 131.77, 123.93, 115.30, 105.97, 98.95, 93.40, 88.01, 83.33, 78.48, 78.24, 77.95, 77.71]
exp_Et = [4.74, 5.21, 6.16, 7.27, 11.14, 18.88, 22.67, 27.73, 31.05, 37.13, 39.58, 41.63, 43.13, 43.14, 43.21, 43.61, 43.87, 46.37]


# Create matrices
Experimental_data = permutedims([ time_data cells./maximum(cells)  exp_TS./maximum(exp_TS) exp_Et./maximum(exp_Et)])
scatter(time_data, cells, label="Cells", xlabel="Time (h)", ylabel="Cells (g/L)", title="Cells vs Time")
scatter(time_data, exp_TS, label="TS", xlabel="Time (h)", ylabel="TS (g/L)", title="TS vs Time")
scatter(time_data, exp_Et, label="EtOH", xlabel="Time (h)", ylabel="EtOH (g/L)", title="EtOH vs Time")


# analysis of the specific growth rate


a = specific_gr_evaluation(permutedims([ time_data cells]),3)
b = specific_gr_evaluation(permutedims([ time_data exp_TS]),3)
c = specific_gr_evaluation(permutedims([ time_data exp_Et]),3)
scatter(a,xlabel="Time (h)",ylabel="Specific growth",label="cells")
plot!(a,xlabel="Time (h)",ylabel="Specific growth",label="cells" )
scatter!(b,xlabel="Time (h)",ylabel="Specific growth",label="Substrate")
plot!(b,xlabel="Time (h)",ylabel="Specific growth",label="Substrate")
scatter!(c,xlabel="Time (h)",ylabel="Specific growth",label="EtOH")
plot!(c,xlabel="Time (h)",ylabel="Specific growth",label="EtOH")
hline!([0.0])
#correlation between the specific growth rate of the cells and the substrate

scatter(a,b,title=string(cor(a,b)),xlabel="Specific growth rate of the cells",ylabel="Specific growth rate of the substrate")
scatter(a,c,title=string(cor(a,c)),xlabel="Specific growth rate of the cells",ylabel="Specific growth rate of the EtOH")  
scatter(b,c,title=string(cor(b,c)),xlabel="Specific growth rate of the substrate",ylabel="Specific growth rate of the EtOH")
#generating the model to fit 
# using inhibition from "Generalization of monod kinetics for analysis of growth data with substrate inhibition"
function Monod_xx_batch(du, u, param, t)
    B1, S, EtOH = u
    K_G,K_E, mu_max_1,Y_1,Y_2,lag= param


  
    # Monod equation
    Q = mu_max_1 *(S)/(K_G + S)
    inhibition =(1 - EtOH/K_E)
    lag_factor = t/(t+lag)

    mu = Q * inhibition * lag_factor

    du[1] =   mu *  B1    #  bio mass dynamics

    du[2] =   -  mu *  B1* Y_1     # Substrate dynamics

    du[3] =    mu * B1 * Y_2   #  EtOH dynamics

  
end





u0 = [Experimental_data[2,1],Experimental_data[3,1],Experimental_data[4,1]]  
param = [ 0.0010361,  0.1866 ,  0.3, 0.0465529, 0.022205,  23.9381 ] #  
ub_p = [500000.5, 1.0, 200.0,200.0, 500.5, 24.0]
lb_p = [0.0,    00.0,  0.0, 0.0, 0.0, 0.0]

test_sim = ODEs_system_sim(Monod_xx_batch,u0, time_data[1],time_data[end], 0.1,param)
plot(test_sim.t, reduce(hcat,test_sim.u)[1,:], label="Cells", xlabel="Time (h)", ylabel="Cells (g/L)")
plot(test_sim.t, reduce(hcat,test_sim.u)[2,:], label="S", xlabel="Time (h)", ylabel="S (g/L)")
plot(test_sim.t, reduce(hcat,test_sim.u)[3,:], label="EtOh", xlabel="Time (h)", ylabel="EtOH (g/L)")



fit_prod_EtOH= Kinbiont.fit_ODEs_System(
    Experimental_data,
    "Test_EtOH",  # Label for dataset
    Monod_xx_batch,  # Custom ODE function
    param,  # Initial parameter guess
    u0 ; # Initial conditions
    ub=ub_p,
    lb=lb_p,
  # maxiter=100000,
)
# the fitting param
fit_prod_EtOH[2]
fitting = fit_prod_EtOH[3]
solt_to_plot = reduce(hcat,fitting.u)
scatter(time_data, Experimental_data[2,:], label="Cells", xlabel="Time (h)", ylabel="Cells (g/L)", title="Cells vs Time")
plot!(time_data, solt_to_plot[1,:], label="Cells", xlabel="Time (h)", ylabel="Cells (g/L)", title="Cells vs Time")

scatter(time_data,  Experimental_data[3,:], label="S", xlabel="Time (h)", ylabel="S (g/L)", title="Substrate vs Time")
plot!(time_data, solt_to_plot[2,:], label="S", xlabel="Time (h)", ylabel="S (g/L)", title="Substrate vs Time")

scatter(time_data,  Experimental_data[4,:], label="EtOH", xlabel="Time (h)", ylabel="EtOH (g/L)", title="EtOH vs Time")
plot!(time_data, solt_to_plot[3,:], label="EtOH", xlabel="Time (h)", ylabel="EtOH (g/L)", title="EtOH vs Time")

