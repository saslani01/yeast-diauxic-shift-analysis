from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import numpy as np

def cluster(gene_expression, min_clusters, max_clusters, auto_stop, thresh):
    best_result = None
    best_silhouette = -1
    prev_silhouette = None
    for n in range(min_clusters, max_clusters + 1):
        k_means = KMeans(n_clusters=n, random_state=13, n_init=10) 
        labels = k_means.fit_predict(gene_expression)
        silhouette = silhouette_score(gene_expression, labels)

        if silhouette > best_silhouette:
            best_result = (n, silhouette, labels)
            best_silhouette = silhouette
        
        print(f"{n}-cluster: silhouette_score={silhouette}")

        if auto_stop and prev_silhouette is not None:
            d = abs(silhouette - prev_silhouette)
            if d < thresh:
                print(f"\nAuto-Stopped at {n} cluster where delta={d} < thresh={thresh}") 
                break
        prev_silhouette = silhouette
    
    print("Number of Cluster =", best_result[0])
    print(f"Silhouette of {best_result[0]} Clusters =", best_silhouette)
    print()

    return best_result

def std_metric(gene_expression):
    return np.std(gene_expression)

def proposed_metric_1(gene_expression):
    before_shift = gene_expression[:4]
    after_shift = gene_expression[4:]
    
    before_shift_mean = np.mean(before_shift)
    after_shift_mean = np.mean(after_shift)
    
    return np.abs(before_shift_mean - after_shift_mean)

def proposed_metric_2(gene_expression):
    return np.max(np.abs(gene_expression[1:] - gene_expression[0]))

def range_metric(gene_expression): 
    return np.max(gene_expression) - np.min(gene_expression)

def top_variable_genes_indicies(gene_expression, top, metric_func):
    
    scores = np.apply_along_axis(metric_func, axis=1, arr=gene_expression)
    return np.argsort(scores)[::-1][:top]
 
def jaccard_coeff(set1, set2):
    intersects = set1 & set2
    union = set1 | set2
    return len(intersects) / len(union), intersects
