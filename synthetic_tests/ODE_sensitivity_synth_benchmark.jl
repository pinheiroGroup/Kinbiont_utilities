
using Kinbiont
using Plots
using Distributions
using Statistics
n_step_sensitivity = 5

# starting time
min_t = 0.0
#
max_t = 250.0
delta_t = 1.0
# size or std of noise
#noise_levels = [0.0001, 0.01, 0.015, 0.020, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05]
model = "aHPM"
noise = 0.05
ode_lb = [0.01, 0.0, 0.3, 0.1]
ode_ub = [0.1, 0.1, 2.5, 10.0]

u0_sim = [0.1]
# matrix of results



u0_sim = [0.1, 0.0]
param_gt = [0.08, 0.02, 1.0, 1.1]
# simulation of data
test_sim = ODE_sim(model, #string of the model
    u0_sim, # starting condition
    min_t, # start time of the sim
    max_t, # final time of the sim
    delta_t, # delta t for poisson approx
    param_gt; # parameters of the ODE model
)
time_sol = reduce(hcat, test_sim.t)
sol_fin = reduce(hcat, test_sim.u)
sol_fin = sum(sol_fin, dims=1)

data_sim = Matrix(vcat(time_sol, sol_fin))

plot(data_sim[1, :], data_sim[2, :])


# adding noise & removing negative values if presents
noise_unifom = rand(Normal(0, noise / 2), length(time_sol))

data_sim[2, :] = data_sim[2, :] .+ noise_unifom


data_sim, index_not_zero = remove_negative_value(data_sim[2, :])

data_sim = Matrix(transpose(hcat(time_sol[index_not_zero], data_sim)))
# fit ODE
sensitivity_test = one_well_morris_sensitivity(data_sim, # dataset first row times second row OD
    "test", # name of the well
    "test_sensitivity", #label of the experiment
    model, # ode model to use 
    ode_lb, # lower bound param
    ode_ub; # upper bound param
    N_step_morris=n_step_sensitivity,
)

