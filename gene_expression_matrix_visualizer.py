import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from analyzer import cluster
import matplotlib.colors as mcolors


def create_heatmap(gene_expression, filename, clustered=False, min_clusters=None, max_clusters=None, clustering_thresh=None):

    if clustered:
        _, _, best_labels = cluster(gene_expression, min_clusters, max_clusters, True, clustering_thresh)

        sorted_indicies = np.argsort(best_labels)
        data = gene_expression[sorted_indicies]
        
        title = "Clustered Gene Expression Matrix"
        ylabel = "Clustered Genes"
    
    else:
        data = gene_expression
        title = "Unclustered Gene Expression Matrix"
        ylabel = "Genes"

    plt.figure(figsize=(8, 10))
    sns.heatmap(
        data,
        cmap = mcolors.LinearSegmentedColormap.from_list("GnBkRd", ["green", "black", "red"]),
        yticklabels=False,
        center=0,
        robust=True
    )
    plt.xlabel("Time Points")
    plt.ylabel(ylabel)
    plt.title(title)

    os.makedirs("outputs", exist_ok=True)
    plt.savefig(f"outputs/{filename}.png")
    plt.close()
