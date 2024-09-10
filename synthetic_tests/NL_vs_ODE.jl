
using Kinbiont
using Plots
using Distributions
using Statistics

n_of_test = 100;

# starting time
min_t = 0.0
#
max_t = 250.0
delta_t = 1.0
# size or std of noise
noise_levels = [0.0001, 0.01, 0.015, 0.020, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05]
model = "logistic"

ode_lb = [0.01, 0.3]
ode_ub = [0.1, 2.5]
ode_u0 = [0.05, 1.0]

NL_lb = [0.3, 0.01, -5 * max_t]
NL_ub = [2.5, 0.1, 5 * max_t]
NL_u0 = [1.0, 0.05, 0.0]


u0_sim = [0.1]
# matrix of results




res_GT = ["rep", "p1", "p2"]

res_ODE = initialize_df_results(model)
res_ODE = vcat(res_ODE, res_GT)
res_ODE = vcat(res_ODE, "noise")

res_NL = initialize_df_results_ode_custom(NL_lb)
res_NL = vcat(res_NL, res_GT)
res_NL = vcat(res_NL, "noise")


# generative model is an ODE

for rep in 1:n_of_test
    println(string("rep =  ", rep))

    # getting random parameters

    param_gt = [rand(Uniform(ode_lb[i], ode_ub[i]), 1) for i in eachindex(ode_ub)]
    param_gt = reduce(hcat, param_gt)[1:end]
    u0_sim = rand(Uniform(0.1, param_gt[2] / 2), 1)
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

    #storing gt
    res_GT_temp = param_gt
    res_GT_temp = vcat(rep, res_GT_temp)

    for noise in noise_levels
        # adding noise & removing negative values if presents
        noise_unifom = rand(Normal(0, noise / 2), length(time_sol))

        data_sim[2, :] = data_sim[2, :] .+ noise_unifom


        data_sim, index_not_zero = remove_negative_value(data_sim[2, :])

        data_sim = Matrix(transpose(hcat(time_sol[index_not_zero], data_sim)))
        # fit ODE
        res_ODE_temp = fitting_one_well_ODE_constrained(
            data_sim, # dataset first row times second row OD
            string(rep), # name of the well
            string("noise_", noise), #label of the experiment
            "logistic", # ode model to `use 
            ode_u0;
            lb=ode_lb, # lower bound param
            ub=ode_ub)
            

        res_ODE_temp = vcat(res_ODE_temp[2], res_GT_temp)
        res_ODE_temp = vcat(res_ODE_temp, noise)
        res_ODE = hcat(res_ODE, res_ODE_temp)
        # fit NL
        res_NL_temp = fit_NL_model(data_sim, # dataset first row times second row OD
            string(rep), # name of the well
            string("noise_", noise), #label of the experiment
            "NL_logistic", #  model to use
            NL_u0;
            lb=NL_lb, # lower bound param
            ub=NL_ub # upper bound param
        )
        res_NL_temp = vcat(res_NL_temp[2], res_GT_temp)
        res_NL_temp = vcat(res_NL_temp, noise)


    
    end
end

# Generative model is NL

res_GT = ["rep", "p1", "p2"]

res_ODE = initialize_df_results(model)
res_ODE = vcat(res_ODE, res_GT)
res_ODE = vcat(res_ODE, "noise")

res_NL = initialize_df_results_ode_custom(NL_lb)
res_NL = vcat(res_NL, res_GT)
res_NL = vcat(res_NL, "noise")

# generative model is an NL

for rep in 1:n_of_test
    println(string("rep =  ", rep))

    # getting random parameters

    param_gt = [rand(Uniform(NL_lb[i], NL_ub[i]), 1) for i in eachindex(ode_ub)]
    param_gt = reduce(hcat, param_gt)[1:end]
    param_gt = push!(param_gt, 0.0)

    # simulation of data
    seq_t = min_t:delta_t:max_t
    time_sol = [seq_t[i] for i in eachindex(seq_t)]
    test_sim = NL_model_logistic(param_gt, time_sol)

    data_sim = Matrix(transpose(hcat(time_sol, test_sim)))
    #storing gt
    res_GT_temp = param_gt[1:(end-1)]
    res_GT_temp = vcat(rep, res_GT_temp)

    for noise in noise_levels
        # adding noise & removing negative values if presents
        noise_unifom = rand(Normal(0, noise / 2), length(data_sim[2, :]))

        data_sim[2, :] = data_sim[2, :] .+ noise_unifom


        data_sim, index_not_zero = remove_negative_value(data_sim[2, :])

        data_sim = Matrix(transpose(hcat(time_sol[index_not_zero], data_sim)))
        # fit ODE
        res_ODE_temp =  fitting_one_well_ODE_constrained(data_sim, # dataset first row times second row OD
        string(rep), # name of the well
        string("noise_", noise), #label of the experiment
        "logistic", # ode model to `use 
        ode_u0;
        lb=ode_lb, # lower bound param
        ub=ode_ub)


        res_ODE_temp = vcat(res_ODE_temp[2], res_GT_temp)
        res_ODE_temp = vcat(res_ODE_temp, noise)
        res_ODE = hcat(res_ODE, res_ODE_temp)
        # fit NL
        res_NL_temp = fit_NL_model(data_sim, # dataset first row times second row OD
        string(rep), # name of the well
        string("noise_", noise), #label of the experiment
        "NL_logistic", #  model to use
        NL_u0;
        lb=NL_lb, # lower bound param
        ub=NL_ub # upper bound param
    )
    res_NL
        res_NL_temp = vcat(res_NL_temp[2], res_GT_temp)
        res_NL_temp = vcat(res_NL_temp, noise)
        res_NL = hcat(res_NL, res_NL_temp)
       

    end
end
