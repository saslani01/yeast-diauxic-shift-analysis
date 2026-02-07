from loader import *
from analyzer import *
from gene_expression_matrix_visualizer import *
from pval_intuition import *
import os

def main():
    # gene_expression refers to the numerical part of the matrix used for the analysis
    # finding the most 228 variable genes (not 230):
    #   * YDR258C appears twice in authors work, so we average their ratios and keep only one
    #   * similarly, when loading 6152 genes from diauxic_raw, grouping by ORF and averaging ratios remove 51 duplicates
    #   * therefore, it makes sence to find min(228, 6152) varaible genes
    df_authors_228, gene_expression_authors_228 = load("./data/230_authors.txt", log2=False)
    df_all, gene_expression_all = load("./data/diauxic_raw_ratios.txt", log2=True)
    min_clusters = 2
    max_clusters = 50
    
    N = 6000
    K = 4260
    n = 50 
    observed_val = 38
    t_vals=[1, 10, 100, 1000, 10000, 1000000]
    t = 1000000 # found to be good number of trials using the plot

    clustering_thresh = 0.005
    alpha_for_p_val = 0.05




    ############################################ most variable genes #########################################
    print(f"{'*' * 50}YEAST_GENE_EXPRESSION_ANALYSIS{'*' * 50}\n")
    print("MOST VARIABLE GENES\n-------------------")
    authors_most_variable_orfs = df_authors_228.index
    metrics = [proposed_metric_1, proposed_metric_2, std_metric, range_metric]
    best_metric = None
    max_overlap = 0
    best_variable_orfs = None
    for m in metrics:
        my_most_variable_orfs = df_all.index[top_variable_genes_indicies(gene_expression_all, 228, m)]
        jaccard, overlap = jaccard_coeff(set(my_most_variable_orfs), set(authors_most_variable_orfs))

        if len(overlap) > max_overlap:
            best_metric = m
            max_overlap = len(overlap)
            best_variable_orfs = my_most_variable_orfs

        print(f"metric={m.__name__:<17} | jaccard={jaccard} | overlap={len(overlap)}")
    print()
    # feed the written file to ontolgy website
    assert best_metric is not None, "No metrics ran successfully!"
    os.makedirs("outputs", exist_ok=True)
    with open(f"./outputs/my_most_variable_orfs_using_{best_metric.__name__}.txt", 'w') as f:
        for orf in best_variable_orfs:
            f.write(f"{orf}\n")
    print(f"Most variable genes were extracted and saved at my_most_variable_orfs_using_{best_metric.__name__}.txt :D")
    print()
    ############################################ Gene Expression Heatmanps ##############################################
    print("CREATING GENE EXPRESSION HEATMAPS\n---------------------------------")
    # gene expression heatmaps
    print(f"Clustering 228_authors_selection, testing for range of min_clusters = {min_clusters} and max_clusters = {max_clusters}")
    create_heatmap(gene_expression_authors_228, "228_clustered", True, min_clusters, max_clusters, clustering_thresh)
    create_heatmap(gene_expression_authors_228, "228_unclustered")
    create_heatmap(gene_expression_all, "all_genes_unclustered")
    print(f"Clustering all_genes, testing for range of min_clusters = {min_clusters} and max_clusters = {max_clusters}")
    create_heatmap(gene_expression_all, "all_genes_clustered", True, min_clusters, max_clusters, clustering_thresh)
    #******************************************************************************************************************#
    print(f"{'*' * 50}P_VALUE_INTUITION_SIMULATIONS{'*' * 50}\n")
    # p-value delta
    print("EMPIRICAL P VALUES FOR (1, 10, ..., 1000000) TRIALS\n---------------------------------------------------")
    print("Generating empirical p-value trial count convergence plot")
    print("N=6000, K=4260, n=50, observed_val=38, t_vals=[1, 10, 100, 1000, 10000, 1000000]")
    pvals, pval_changes = create_p_val_trial_count_convergence_plot(N, K, n, observed_val, t_vals)
    print("pvals:", pvals)
    print("pval_changes:", pval_changes)
    print()
    ############################################ P-value Intuition (how GO works) ############################################
    print("P-VALUES CALCULATED USIGN DIFFERENT METHODS\n-------------------------------------------")
    emp_pval, bin_pval, hyp_pval = p_vals_using_different_models(N, K, n, observed_val, t=1000000)
    print(f"Empirical p-val: {emp_pval:.4f}")
    print(f"Binomial p-val: {bin_pval:.4f}")
    print(f"Hypergeom p-val: {hyp_pval:.4f}")
    print()

    # different observed_va;ss produce different pvals (using alpha = 0.05)
    print(f"DIFFERENET observed_val PRODUCED DIFFERENT P VALUES (ALPHA = {alpha_for_p_val})\n----------------------------------------------------------------")
    
    obs_35_pval = empirical_pval(N, K, n, t, observed_val=35)
    obs_38_pval = empirical_pval(N, K, n, t, observed_val=38)
    obs_40_pval = empirical_pval(N, K, n, t, observed_val=40)

    obs_35_statitical_sig = is_significant(obs_35_pval, thresh=alpha_for_p_val)
    obs_38_statitical_sig = is_significant(obs_38_pval, thresh=alpha_for_p_val)
    obs_40_statitical_sig = is_significant(obs_40_pval, thresh=alpha_for_p_val)
    print(f"observed_val of {35} (pval={obs_35_pval}) has significance = {obs_35_statitical_sig}")
    print(f"observed_val of {38} (pval={obs_38_pval}) has significance = {obs_38_statitical_sig}")
    print(f"observed_val of {40} (pval={obs_40_pval}) has significance = {obs_40_statitical_sig}")
    print()

    # finding the k and K such that p_val lie in p_min and p_max
    print("FINDING (SUCCESS_SAMPLE, SUCCESS_POPULATION) SUCH THAT P VALUE WITHIN RANGE 0.003 and 0.01\n------------------------------------------------------------------------------------------")
    success_in_population, success_in_sample, pval = find_sucess_sample_and_population_params_for_pval_range(N, n, p_min=0.003, p_max=0.01)
    print(f"success_in_population={success_in_population} success_in_sample={success_in_sample} pval={pval}")
    
if __name__ == "__main__":
    main()
