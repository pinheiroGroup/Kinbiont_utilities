using Kinbiont
using CSV
using Tables
using StatsBase

################################

ub_exp = [0.3]
lb_exp = [00.00]

ub_logistic = [0.1, 0.9]
lb_logistic = [0.0001, 0.000001]

ub_alogistic = [0.1, 0.9, 10.0]
lb_alogistic = [0.0001, 0.000001, 0.001]

ub_hpm = [0.3, 0.8, 0.9]
lb_hpm = [0.0001, 0.000001, 0.001]


ub_hpm_exp = [0.1, 1.0]
lb_hpm_exp = [0.0001, 0.0000001]



list_of_models = ["HPM", "HPM_exp", "logistic", "exponential"]
list_ub_param = [ub_hpm, ub_hpm_exp, ub_logistic, ub_exp]
list_lb_param = [lb_hpm, lb_hpm_exp, lb_logistic, lb_exp]
list_guess = list_lb_param .+ (list_ub_param .- list_lb_param)./2

noise_uniform = zeros(9)
segment_number = [0,1, 2, 3, 4, 5]
data_1 = CSV.File("/data/diauxic/data_modified_1.csv", transpose=true, missingstring=nothing)



path_to_results = string("/res/diauxic/")
path_to_plot = string("/plot/diauxic/")

# number of segment
segment_number = [0, 1, 2, 3 ,4,5]

x_size=500
y_size =700
guidefontsize=18
tickfontsize=16
legendfontsize=10

names_of_cols = propertynames(data_1)
data_temp = data_1[names_of_cols[1]]
times_data = convert(Vector{Float64}, data_temp) ./ 60


for nseg in segment_number

    # 


    for r in names_of_cols[2:end]

        data_OD = Matrix(transpose(hcat(times_data, data_1[r])))


        param_opt = segmentation_ODE(
            data_OD, # dataset first row times second row OD
            string(r), # name of the well
            "data_1", #label of the experiment
            list_of_models, # ode models to use
            list_guess,
            nseg;
            lb_param_array=list_lb_param, # lower bound param
            ub_param_array=list_ub_param, # upper bound param
            detect_number_cpd=false,
            fixed_cpd=false,
            type_of_detection="slinding_win",
            type_of_curve="deriv",
            pt_avg=3, # number of the point to generate intial condition
            smoothing=false, # the smoothing is done or not?            win_size=10, # numebr of the point to generate intial condition
            pt_smooth_derivative=6,
            win_size = 10
            )



   # writing results
    # writng optim res
    CSV.write(string(path_to_results,"/param_diauxic_",string(r),"_",nseg,".csv"),Tables.table(param_opt[2] ))

    #writing fits
    CSV.write(string(path_to_results,"/fits/fits_diauxic_",string(r),"_",nseg,".csv"),Tables.table( param_opt[3]))
    # writing data pre-processed
    CSV.write(string(path_to_results,"/data/times_diauxic_",string(r),"_",nseg,".csv"),Tables.table(param_opt[4] ))
    # writing AIC
    CSV.write(string(path_to_results,"/AIC/AIC_diauxic_",string(r),"_",nseg,".csv"),Tables.table([nseg string(r) param_opt[end]]) )
    # writing cp
    CSV.write(string(path_to_results,"/cp/CP_diauxic_",string(r),"_",nseg,".csv"),Tables.table(param_opt[5] ))

    # plotting
    display(
        Plots.scatter(
            data_OD[1,:],
            data_OD[2,:],
            xlabel="Time",
            ylabel="Arb. Units",
            label=["Data " nothing],
            markersize=2,
            color=:black,
            title=string(string(r)),
            guidefontsize=guidefontsize,
            tickfontsize=tickfontsize,
            legendfontsize=legendfontsize,
            size=(y_size,x_size),
        ),
        )

        display(
            Plots.plot!(
                param_opt[4],
                param_opt[3],
                xlabel="Time",
                ylabel="Arb. Units",
                label=[string("Fitting ") nothing],
                c=:red,
                guidefontsize=guidefontsize,
                tickfontsize=tickfontsize,
                legendfontsize=legendfontsize,
                size=(y_size,x_size),
            ),
        )


        display( Plots.vline!(
            [param_opt[5]],
            c=:black,
            label=[string("Change point") nothing],
            guidefontsize=guidefontsize,
            tickfontsize=tickfontsize,
            legendfontsize=legendfontsize,
            size=(y_size,x_size),
            ) )


            png(string(path_to_plot,  "segmented_fit_", string(r),"_",nseg, ".png"))



    end


end




 