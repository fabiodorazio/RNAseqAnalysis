import scanpy as sc

def clustering(matrix, n_principal_component):
    '''
    scanpy returns the variability of each gene in the dataset by assigning a dispersion score
    dispersion can be visualised using `sc.pl.highly_variable_genes(matrix)`
    '''
    # find the 2000 most variable genes
    sc.pp.highly_variable_genes(matrix, n_top_genes = 2000)
    # select only the highly variable genes
    matrix = matrix[:, matrix.var.highly_variable]
    # regress out the differences that arise due to total count, mito and ribo
    sc.pp.regress_out(matrix, ['total_counts', 'pct_counts_mt', 'pct_counts_ribo'])
    # normalise each gene to the unit variance of that gene
    sc.pp.scale(matrix, max_value=10)
    # reduce the dimesion of the datasets using PCA. Default is 50 PCs
    sc.tl.pca(matrix, svd_solver='arpack')

    # umap. select number of PCs based on the elbow plot
    sc.pp.neighbors(matrix, n_pcs = n_principal_component)
    # this projects the data from the 30 PC dimensions to 2 dimensions in the Umap
    sc.tl.umap(matrix)

    return()