
# Modeling and kinetic study of bio-ethanol production from soy protein  zconcentrate by-product
using Kinbiont
using Plots
# Data from 
# Time points
time_data = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 28.0, 32.0, 36.0, 40.0, 44.0]

# Experimental data
cells = exp.([7.51, 7.51, 7.59, 7.66, 7.75, 7.89, 7.96, 8.08, 8.13, 8.14, 8.14, 8.15, 8.15, 8.15, 8.16, 8.16, 8.16, 8.16])
exp_TS = [170.20, 165.62, 161.15, 159.55, 150.67, 140.66, 131.77, 123.93, 115.30, 105.97, 98.95, 93.40, 88.01, 83.33, 78.48, 78.24, 77.95, 77.71]
exp_Et = [4.74, 5.21, 6.16, 7.27, 11.14, 18.88, 22.67, 27.73, 31.05, 37.13, 39.58, 41.63, 43.13, 43.14, 43.21, 43.61, 43.87, 46.37]


# Create matrices
Experimental_data = permutedims([ time_data cells  exp_TS exp_Et])
scatter(time_data, cells, label="Cells", xlabel="Time (h)", ylabel="Cells (g/L)", title="Cells vs Time")
scatter(time_data, exp_TS, label="TS", xlabel="Time (h)", ylabel="TS (g/L)", title="TS vs Time")
scatter(time_data, exp_Et, label="EtOH", xlabel="Time (h)", ylabel="EtOH (g/L)", title="EtOH vs Time")



# generating the model to fit 

function Growth_with_EtOH_prod(du, u, param, t)
    B1, G, EtOH = u
    k_1,k_2,K_G,K_E, mu_max_1,mu_max_2 ,D,G_R = param


  
    # Monod equation
    Q = k_1 *(G)/(K_G + G)
    Q_M = k_2 *(EtOH)/(K_E + EtOH)




    du[1] =  mu_max_1 * Q * B1 + Q_M *   mu_max_2 * B1 - D*B1# - EtOH_prod * du[1] # bio mass dynamics

    du[2] =   D * (G-G_R)- Q * B1 - D*B1   # Substrate dynamics

    du[3] = mu_max_1 * Q * B1 - mu_max_2 *Q_M * B1 - D*B1   #  EtOH dynamics

  
end





u0 = [Experimental_data[2,1],Experimental_data[3,1],Experimental_data[4,1]]  # Initial conditions NOTE THE DOUBLE OF THE INITIAL GLUCOSE
param = [2.5, 0.03 ,0.02, 0.1,2.5, 0.03 ,0.02, 0.1]  #  
ub_p = [500.5, 300000.0,96.0,10.0,500.5, 300000.0,96.0,10.0]
lb_p = [0.0, 0.0,0.0,0.0,0.0, 0.0,0.0,0.0]


fit_prod_EtOH= Kinbiont.fit_ODEs_System(
    Experimental_data,
    "Test_EtOH",  # Label for dataset
    Growth_with_EtOH_prod,  # Custom ODE function
    param,  # Initial parameter guess
    u0 ; # Initial conditions
    ub=ub_p,
    lb=lb_p,
  #  maxiter=100,
)

fitting = fit_prod_EtOH[3]
solt_to_plot = reduce(hcat,fitting.u)
scatter(time_data, cells, label="Cells", xlabel="Time (h)", ylabel="Cells (g/L)", title="Cells vs Time")
plot!(time_data, solt_to_plot[1,:], label="Cells", xlabel="Time (h)", ylabel="Cells (g/L)", title="Cells vs Time")

scatter(time_data, exp_TS, label="TS", xlabel="Time (h)", ylabel="TS (g/L)", title="TS vs Time")
plot!(time_data, solt_to_plot[2,:], label="Cells", xlabel="Time (h)", ylabel="Cells (g/L)", title="Cells vs Time")

scatter(time_data, exp_Et, label="EtOH", xlabel="Time (h)", ylabel="EtOH (g/L)", title="EtOH vs Time")
plot!(time_data, solt_to_plot[3,:], label="Cells", xlabel="Time (h)", ylabel="Cells (g/L)", title="Cells vs Time")

