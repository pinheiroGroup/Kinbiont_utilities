using Kinbiont
using Plots
using CSV, DataFrames
using Statistics
using DelimitedFiles
using Random




Kinbiont_res_test = readdlm("/Users/fabrizio.angaroni/Documents/KinBiont_utilities-main/Fig_4_5/df_for_ML/res_clean_ML_richards.csv", ',')
annotation_test = readdlm("/Users/fabrizio.angaroni/Documents/KinBiont_utilities-main/Fig_4_5/df_for_ML/annotation_clean_richards.csv", ',')

ordered_strain = annotation_test[:, end]
n_folds = 10
feature_names = unique(annotation_test[1, 2:end])[2:(end-1)]

seq_of_depths = [ -1,2,3,4,5,6,7,8]


gr_n_folds_cv = vcat("Strain","depth", [i for i in 1:n_folds])
avg_gr_n_folds_cv = vcat("Strain","depth", "avg_R2")
gr_impurity_importance = vcat("Strain","depth", feature_names)
gr_split_importance = vcat("Strain","depth", feature_names)




nmax_n_folds_cv = vcat("Strain","depth", [i for i in 1:n_folds])
avg_nmax_n_folds_cv = vcat("Strain","depth", "avg_R2")
nmax_impurity_importance = vcat("Strain","depth", feature_names)
nmax_split_importance = vcat("Strain","depth", feature_names)


lag_n_folds_cv = vcat("Strain","depth", [i for i in 1:n_folds])
avg_lag_n_folds_cv = vcat("Strain","depth", "avg_R2")
lag_impurity_importance = vcat("Strain","depth", feature_names)
lag_split_importance = vcat("Strain","depth", feature_names)

list_strain = unique(Kinbiont_res_test[(3), 2:end])



seed = Random.seed!(12384)


for d in seq_of_depths

    for s in list_strain

         println("#############################")
         println(s)

        index_strain = findall(s .== ordered_strain)
        feature_matrix = annotation_test[index_strain, 2:(end-1)]
        Kinbiont_results = Kinbiont_res_test[:, index_strain]
        dt_gr = downstream_decision_tree_regression(Kinbiont_results,
            feature_matrix,
            9;
            do_pruning=false,
            pruning_accuracy=1.00,
            verbose=true,
            do_cross_validation=true,
            max_depth=d,
            n_folds_cv=n_folds,
            seed=seed
        )


        dt_nmax = downstream_decision_tree_regression(Kinbiont_results,
            feature_matrix,
            5;
            do_pruning=false,
            pruning_accuracy=1.00,
            verbose=true,
            do_cross_validation=true,
            max_depth=d,
            n_folds_cv=n_folds,
            seed=seed)


        dt_lag = downstream_decision_tree_regression(Kinbiont_results,
            feature_matrix,
            8;
            do_pruning=false,
            pruning_accuracy=1.00,
            verbose=true,
            do_cross_validation=true,
            max_depth=d,
            n_folds_cv=n_folds,
            seed=seed)

        temp_gr_n_folds_cv = vcat(d, dt_gr[4])
        temp_gr_n_folds_cv = vcat(s, temp_gr_n_folds_cv)
        gr_n_folds_cv = hcat(gr_n_folds_cv, temp_gr_n_folds_cv)

        temp_gr_impurity_importance = vcat(d, dt_gr[2])
        temp_gr_impurity_importance = vcat(s, temp_gr_impurity_importance)
        gr_impurity_importance = hcat(gr_impurity_importance, temp_gr_impurity_importance)

        temp_gr_split_importance = vcat(d, dt_gr[3])
        temp_gr_split_importance = vcat(s, temp_gr_split_importance)
        gr_split_importance = hcat(gr_split_importance,temp_gr_split_importance)

        temp_nmax_n_folds_cv= vcat(d, dt_nmax[4])
        temp_nmax_n_folds_cv = vcat(s, temp_nmax_n_folds_cv)
        nmax_n_folds_cv = hcat(nmax_n_folds_cv,temp_nmax_n_folds_cv)

        temp_nmax_impurity_importance= vcat(d, dt_nmax[2])
        temp_nmax_impurity_importance = vcat(s, temp_nmax_impurity_importance)
        nmax_impurity_importance = hcat(nmax_impurity_importance, temp_nmax_impurity_importance)

        temp_nmax_split_importance= vcat(d, dt_nmax[3])
        temp_nmax_split_importance = vcat(s, temp_nmax_split_importance)
        nmax_split_importance = hcat(nmax_split_importance, temp_nmax_split_importance)

        temp_lag_n_folds_cv = vcat(d, dt_lag[4])
        temp_lag_n_folds_cv = vcat(s, temp_lag_n_folds_cv)
        lag_n_folds_cv = hcat(lag_n_folds_cv,temp_lag_n_folds_cv)

        temp_lag_impurity_importance = vcat(d, dt_lag[2])
        temp_lag_impurity_importance = vcat(s, temp_lag_impurity_importance)
        lag_impurity_importance = hcat(lag_impurity_importance, temp_lag_impurity_importance)

        temp_lag_split_importance = vcat(d, dt_lag[3])
        temp_lag_split_importance = vcat(s, temp_lag_split_importance)
        lag_split_importance = hcat(lag_split_importance, temp_lag_split_importance)

    end

end



CSV.write("gr_split_imp_clean_loss.csv", Tables.table(gr_split_importance))
CSV.write("nmax_split_imp_clean_loss.csv", Tables.table(nmax_split_importance))
CSV.write("lag_split_imp_clean_loss.csv", Tables.table(lag_split_importance))


CSV.write("gr_impurity_imp_clean_loss.csv", Tables.table(gr_impurity_importance))
CSV.write("nmax_impurity_imp_clean_loss.csv", Tables.table(nmax_impurity_importance))
CSV.write("lag_impurity_imp_clean_loss.csv", Tables.table(lag_impurity_importance))




CSV.write("gr_coeff_clean_loss.csv", Tables.table(gr_n_folds_cv))
CSV.write("nmax_coeff_clean_loss.csv",Tables.table(nmax_n_folds_cv))
CSV.write("lag_coeff_clean_loss.csv", Tables.table(lag_n_folds_cv))