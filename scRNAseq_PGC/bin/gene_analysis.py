import scanpy as sc
import scvi
import pandas as pd
import numpy as np

def cell_frequency(matrix):
    '''
    label each cell based on the sample of origin
    calculates frequency for each sample
    '''
    # get the total number of cells per sample
    num_tot_cells = matrix.obs.groupby(['Sample']).count()
    # creates a dictionary to zip the row names as keys and the total number of cells as values
    num_tot_cells = dict(zip(num_tot_cells.index, num_tot_cells.doublet))
    
    # gets the total number of the group
    cell_type_counts = matrix.obs.groupby(['Sample', 'condition', 'cell_type']).count()
    cell_type_counts = cell_type_counts[cell_type_counts.sum(axis = 1) > 0].reset_index()
    cell_type_counts = cell_type_counts[cell_type_counts.columns[0:4]]
    # use the dictionary to create a new column
    cell_type_counts['total_cells'] = cell_type_counts.Sample.map(num_tot_cells).astype(int)
    # obtain the frequency
    cell_type_counts['frequency'] = cell_type_counts.doublet / cell_type_counts.total_cells
    
    return(matrix)


def diff_expression(matrix):
    '''
    Differential gene expression is modelled using a negative binomial distribution. This is important to account for skewed read counts and batch effects
    '''
    # train the model
    model  = scvi.model.SCVI.load('model.model', matrix)

    scvi_de = model.differential_expression(
    idx1 = [matrix.obs['cell_type'] == 'AT1'],
    idx2 = [matrix.obs['cell_type'] == 'AT2']
    )
    scvi_de = scvi_de[(scvi_de['is_de_fdr_0.05']) & (abs(scvi_de.lfc_mean) > .5)]
    scvi_de = scvi_de.sort_values('lfc_mean')
    scvi_de = scvi_de[(scvi_de.raw_normalized_mean1 > .5) | (scvi_de.raw_normalized_mean2 > .5)]

    return(scvi_de)



