import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc

def plot_doublets(matrix, output_plot_dir):
    sns.displot(matrix[matrix.prediction == 'doublet'], x = 'dif')
    plt.savefig(f'{output_plot_dir}/Figure1.png')
    plt.close()


    
def plot_qc_metrics(matrix, output_plot_dir):
    '''
    Use to get rid of outliers
    Significant higher genes or counts means it is an artifact
    High mitochondrial representation means stressed and apoptotic cells. Set a filter between 5 and 20%
    '''
    sc.pl.violin(matrix, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'], 
             jitter=0.4, multi_panel=True)
    plt.savefig(f'{output_plot_dir}/Figure2.png')
    plt.close()
