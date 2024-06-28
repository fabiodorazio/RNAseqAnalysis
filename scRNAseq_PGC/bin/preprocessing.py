import scanpy as sc
import scvi
import pandas as pd
import numpy as np

def load_matrix_file(file):
    '''
    Load matrix file made of gene ensembls in the rows and individual cells as columns
    '''
    # read matrix and transpose it to align with scanpy requirements: rows as individual cells
    input_matrix = sc.read_csv(file).T
    
    return(input_matrix)


def label_doublets(matrix):
    '''
    removing doublets
    '''
    # filter only genes with expression in at least 10 cells
    sc.pp.filter_genes(matrix, min_cells=10)
    # keep the top highly variable genes that describe the model at best
    sc.pp.highly_variable_genes(matrix, n_top_genes=2000, subset=True, flavor='seurat_v3')
    # train model to predict the doublets
    scvi.model.SCV.setup_anndata(matrix) # model setup with default parameters
    v = scvi.model.SCVI(matrix)
    v.train()
     # train the SOLO model to predict doublets
    solo = scvi.external.SOLO.from_scvi_model(v)
    solo.train()
    # predictions on whether a cell is a doublet or not. The higher the score the more likely the cell is a doublet 
    predictions = solo.predict()
    # make a new column as predicted label
    predictions['prediction'] = solo.predict(soft=False)
    # remove all the -0 at the end of the barcode that scvi adds
    predictions.index = predictions.index.map(lambda x:x[-2])

    # count doublets and singlets. Expected values are between 5 and 20% of doublets
    singlets = predictions.groupby('prediction').count()
    print(f"number of singles it {singlets.shape}")
    # calculate the difference between doublets and singlets to be less stringent on the filter for doublets
    predictions['difference'] = predictions.doublet - predictions.singlet
    # predict doublets for cells that have label doublets and a difference higher than 1
    doublets = predictions[(predictions.prediction == 'doublet') & (predictions.difference > 1)]

    return(doublets)


def filter_out_doublets(matrix, doublets):
    # add doublet column with true or false values based on the match with the doublets dataset
    matrix.obs['doublet'] = matrix.obs.index.isin(doublets.index)
    # filter out predicted doublets from the original table
    singlet_matrix = matrix[~matrix.ob.doublet]
    
    return(singlet_matrix)


def label_mithocondrial_genes(matrix):
    '''
    mitochondrial genes are labelled using MT- at the start
    '''
    matrix.var['mt'] = matrix.var[matrix.var.index.str.startswith('MT-')]

    return(matrix)

def label_ribosomal_genes(matrix, ribo_url): # OPTIONAL
    '''
    import list of ribosomal genes from broad institute
    filter out genes
    '''
    ribo_genes = pd.read_table(ribo_url, skiprows=2, header=None)
    # add column for ribosomal genes on original matrix 
    matrix.var['ribo'] = matrix.var_names.isin(ribo_genes[0].values)

    return(matrix)


def calculate_qc_metrics(matrix, ribo_url):
    '''
    use scanpy to calculate qc metrics
    '''
    # label mitochondrial and ribosomal RNAs
    mito_matrix = label_mithocondrial_genes(matrix)
    ribo_mito_matrix = label_ribosomal_genes(mito_matrix, ribo_url)

    # the new columns added mt and ribo are fed to the function to calculate the qc
    sc.pp.calculate_qc_metrics(ribo_mito_matrix, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
    # sort increasingly by the number of cells per gene
    ribo_mito_matrix.var.sort_values('n_cells_by_counts')
    # remove genes that are not in at least 4 cells
    sc.pp.filter_genes(ribo_mito_matrix, min_cells=4)

    # filter by total counts. Arbitrary based on the data
    sc.pp.filter_cells(ribo_mito_matrix, min_genes=200)

    return(ribo_mito_matrix)



def filtering(matrix, percentile_choice, mito_filter, ribo_filter):
    '''
    filter out cells based on the number of genes by counts
    mitocondrial and ribosomal rna content
    '''
    # calculate upper limit with the 98th percentile. You can set an arbitrary number of genes. Ex = 3000
    U_limit = np.quantile(matrix.obs.n_genes_by_counts.values, percentile_choice)
    # filter genes by U_limit
    matrix = matrix[matrix.obs.n_genes_by_counts < U_limit]

    # filter out cells that have mitochondrial rna reads > 20%
    matrix = matrix[matrix.obs.pct_counts_mt < mito_filter]
    # filter out cells that have ribosomal rna reads > 2
    matrix = matrix[matrix.obs.pct_counts_ribo < ribo_filter]
    
    return(matrix)


def normalisation(matrix):
    '''
    sc sequencing the variation is high between cells
    Normalisation is required for comparison between cells
    There is often an order of magnitude difference between the counts per cell in the dataset
    '''
    # normalize every cell to 10,000 UMI so that the total number fo counts per cell is equal
    sc.pp.normalize_total(matrix, target_sum=1e4)
    # change to log counts
    sc.pp.log1p(matrix)
    # save the new matrix in the raw table in matrix. This populates the raw table
    matrix.raw = matrix
    return(matrix)

