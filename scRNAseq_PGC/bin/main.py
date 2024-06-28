import os
import scanpy as sc
import scvi
import pandas as pd
import numpy as np
import scipy
import argparse
import time

import preprocessing as prp
import clustering as clst
import plots
import utils


# ribosomal RNA list
ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"

# set argument parameters
def get_args():
    parser = argparse.ArgumentParser(prog="scRNAsequencer")
    parser.add_argument("input_sc", 
                        help = "Input sc table")
    parser.add_argument("-n", dest="name", default=None, 
                        help="Output file name without the format: can be the name of the analysis, individual, date, etc...")
    parser.add_argument("-o", dest="output_dir", default="../Outputs", 
                        help="Output directory path")
    parser.add_argument("--filter-percent", dest="pct", default=None, 
                        help="Set the percentile used to filter cells on number of genes")
    parser.add_argument("--filter-ribo", dest="ribo", default=2, type=int, 
                        help="Arbitrary percent value used for filtering out cells based on ribosomal RNA")
    parser.add_argument("--filter-mito", dest="mito", default=20, type=int, 
                        help="Arbitrary percent value used for filtering out cells based on mitochondrial RNA")
    parser.add_argument("--summary-plots", dest="plots", action="store_true", 
                        help="True or false: generate summary plots")
    parser.add_argument("--principal-components", dest="pcs", default=30, type=int,
                        help="Number of Principal Components selected to perform the dimensionality reduction")

    args = parser.parse_args()

    input = args.input_sc
    name = args.name
    output = args.output_dir
    pct = args.pct
    ribo = args.ribo
    mito = args.mito
    plots = args.plots
    pcs = args.pcs

    return (input, name, output, pct, ribo, mito, plots, pcs)


if __name__ == "__main__":
    INPUT, NAME, OUTPUT, PCT, RIBO, MITO, PLOTS, PCS = get_args()
    print("Starting the aligner...")
    time.sleep(2)

    if NAME is None:
        NAME = utils.get_basename(INPUT)

     # checks if output exists or create one
    output_dir = utils.check_output_dir(OUTPUT)
    # creates output path for the aligner
    output_analysis = output_dir + NAME + ".csv"
    # sets location for plot output, checks if exists or creates one
    output_plot_dir = utils.check_output_dir(output_dir + '/Plot_outputs/')
    output_plot_dir_basename = output_plot_dir + NAME


################### RUN FUNCTIONS ###################
    
m_matrix = prp.load_matrix_file("../assets/zebrafish_PGC_rep1.csv")

# remove doublets
predicted_doublets = prp.label_doublets(m_matrix)
singlets = prp.filter_out_doublets(m_matrix, predicted_doublets)

if PLOTS:
    plots.plot_doublets(predicted_doublets, output_plot_dir_basename)

# QC metrics
qc_metrics_table = prp.calculate_qc_metrics(singlets)

if PLOTS:
    plots.plot_qc_metrics(qc_metrics_table, output_plot_dir_basename)

filtered_table = prp.filtering(qc_metrics_table, PCT, MITO, RIBO)

# clustering
clustered = clst.clustering(filtered_table, PCS)

if PLOTS:
    plots.plot_umap(clustered, output_plot_dir_basename)