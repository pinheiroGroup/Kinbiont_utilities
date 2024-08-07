using KinBiont
using Plots
using CSV, DataFrames
using Statistics
using DelimitedFiles
using Random
using DecisionTree
using AbstractTrees
using MLJDecisionTreeInterface
using TreeRecipe

kimchi_res_test = readdlm("/df_for_ML/res_clean_ML_richards.csv", ',')
annotation_test = readdlm("df_for_ML/annotation_clean_richards.csv", ',')

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

list_strain = ["N. soli", "mixture"]


seed = Random.seed!(1234)



for s in list_strain

         println("#############################")
         println(s)

        index_strain = findall(s .== ordered_strain)
        feature_matrix = annotation_test[index_strain, 2:(end-1)]
        kimchi_results = kimchi_res_test[:, index_strain]

        dt_gr = downstream_decision_tree_regression(kimchi_results,
            feature_matrix,
            9;
            do_pruning=false,
            pruning_accuracy=1.00,
            verbose=true,
            do_cross_validation=true,
            max_depth=2,
            n_folds_cv=n_folds,
            seed=seed
        )
        wt = DecisionTree.wrap(dt_gr[1], (featurenames = feature_names,))
 
 
        p2 = Plots.plot(wt, 0.9, 0.2; size = (900,400), connect_labels = ["yes", "no"])
 
        savefig(p2,string(string(s),"_dt_gr.svg"))

        dt_nmax = downstream_decision_tree_regression(kimchi_results,
            feature_matrix,
            5;
            do_pruning=false,
            pruning_accuracy=1.00,
            verbose=true,
            do_cross_validation=true,
            max_depth=2,
            n_folds_cv=n_folds,
            seed=seed)

            wt = DecisionTree.wrap(dt_nmax[1], (featurenames = feature_names,))
 
 
            p2 = Plots.plot(wt, 0.9, 0.2; size = (900,400), connect_labels = ["yes", "no"])
     
            savefig(p2,string(string(s),"_dt_N_max.svg"))
    
        dt_lag = downstream_decision_tree_regression(kimchi_results,
            feature_matrix,
            8;
            do_pruning=false,
            pruning_accuracy=1.00,
            verbose=true,
            do_cross_validation=true,
            max_depth=2,
            n_folds_cv=n_folds,
            seed=seed)

            wt = DecisionTree.wrap(dt_lag[1], (featurenames = feature_names,))
 
 
            p2 = Plots.plot(wt, 0.9, 0.2; size = (900,400), connect_labels = ["yes", "no"])
            savefig(p2,string(string(s),"_dt_lag.svg"))


            #saving leaf distrib

            CSV.write(string(string(s),"_leaves_dt_lag.csv"),Tables.table(permutedims(dt_lag[end])))
            CSV.write(string(string(s),"_leaves_dt_gr.csv"),Tables.table(permutedims(dt_gr[end])))
            CSV.write(string(string(s),"_leaves_dt_nmax.csv"),Tables.table(permutedims(dt_nmax[end])))

end








