import pandas as pd
import numpy as np

def load(path, log2):
    df = pd.read_csv(path, sep='\t', index_col=0) # drop the row number
    ratio_cols = [col for col in df.columns if "Ratio" in col]
    
    # averaging the ratios for rows having the same ORF (grouping by ORF and dropping name aliases)
    # drop rows that dont have ORF id (we actually dont have any na here, but they would useless especially for gene ontology)
    df = df.groupby("ORF", dropna=True)[ratio_cols].mean()
    gene_expression = df.values
    return (df, np.log2(gene_expression)) if log2 else (df, gene_expression)
