using Kinbiont
using CSV
using Statistics
using Distributions
using Plots

n_step_sensitivity = 5

# starting time
min_t = 0.0
#
max_t = 250.0
delta_t = 1.0
# size or std of noise
#noise_levels = [0.0001, 0.01, 0.015, 0.020, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05]
model = ["NL_Richards"]
noise = 0.05
NL_lb = [[0.01, 0.0, 0.01, 0.1]]
NL_ub = [[2.1, 10.1, 0.1, 1000.0]]
param_gt = [[1.0, 1.02, 0.05, 200.1]]





# simulation of data
seq= min_t:delta_t:max_t
time_sol = [seq[i] for i in eachindex(seq)]
test_sim = NL_model_Richards(param_gt[1],time_sol)


data_sim = Matrix(transpose(hcat(time_sol, test_sim)))

plot(data_sim[1, :], data_sim[2, :])


# adding noise & removing negative values if presents
noise_unifom = rand(Normal(0, noise / 2), length(time_sol))

data_sim[2, :] = data_sim[2, :] .+ noise_unifom


data_sim, index_not_zero = remove_negative_value(data_sim[2, :])

data_sim = Matrix(transpose(hcat(time_sol[index_not_zero], data_sim)))

nl_fit =  NL_model_selection(data_sim, # dataset first row times second row OD
"test", 
"test_model_selection",
model, #  model to use
param_gt;
nrep =n_step_sensitivity,
method_of_fitting ="Morris_sensitivity",
lb_param_array = NL_lb,
ub_param_array = NL_ub
)

