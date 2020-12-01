import pandas as pd
import fastcluster
import numpy as np
from Bio import Phylo

def cluster_ks(df):
    """
    Cluster a data frame with Ks estimates using average-linkage hierarchical
    clustering.
    """
    X = pairwise_matrix(df, "dS")
    Y = alc_cluster(X)
    T = to_phylo(Y, X)
    return T

def alc_cluster(df, nanpolicy=1000):
    # fill NaN values with something larger than all the rest, not a
    # foolproof approach, but should be reasonable in most cases
    if np.any(np.isnan(df)):
        logging.warning("Data contains NaN values, replaced by "+str(nanpolicy))
        df.fillna(nanpolicy, inplace=True)
    return fastcluster.average(df)

def pairwise_matrix(df, col):
    d = {}
    for p in df.index:
        g1 = df.loc[p].gene1
        g2 = df.loc[p].gene2
        x  = df.loc[p][col]
        if not g1 in d:
            d[g1] = {}
        if not g2 in d:
            d[g2] = {}
        d[g1][g2] = x
        d[g2][g1] = x
    df = pd.DataFrame.from_dict(d).fillna(0.)
    # df[df.index] ensures the result is a proper (symmetric) distance matrix
    return df[df.index]

def to_phylo(clustering, df):
    n = len(df.index)
    heights = {i: 0 for i in range(n)}
    clades = [Phylo.Newick.Clade(0, df.index[i]) for i in range(n)]
    for row in clustering:
        i = int(row[0])
        j = int(row[1])
        c1 = clades[i]
        c2 = clades[j]
        height = row[2] 
        newnode = Phylo.Newick.Clade(0, n)
        c1.branch_length = height - heights[i]
        c2.branch_length = height - heights[j]
        newnode.clades.extend([c1, c2])
        heights[n] = height
        clades.append(newnode)
        n += 1
    return clades[-1]
    
