import os
import argparse
import pandas as pd
from CellAlignDTW import CellAlignDTW, gam_smooth_expression
from CellAlignDTW.pp import create_cellrank_probability_df

def main():
    parser = argparse.ArgumentParser(description="Run CellAlignDTW analysis")
    parser.add_argument("--path", required=True, help="Path to the directory containing the h5ad files")
    parser.add_argument("--clusters", required=True, help="Comma-separated list of clusters to analyze")
    parser.add_argument("--gene_list", required=True, help="Path to the file containing the list of DE genes")
    args = parser.parse_args()

    print(f"Starting CellAlignDTW analysis with path: {args.path}, clusters: {args.clusters}, and gene list: {args.gene_list}")

    data = """
RU151,NSCLC_1,SCLC-N_2,SCLC-A,NonNE SCLC_1
RU263,NSCLC_2,SCLC-A_3,SCLC-N_1,SCLC-N_3
RU581,NSCLC,SCLC-A_4,SCLC-A_5,SCLC-N
RU831,NSCLC,SCLC-A_5,SCLC-A_6,SCLC-A_3,SCLC-N
RU942,NSCLC_1,NSCLC_3,SCLC-A_3,SCLC-N_2,SCLC-N_1
RU1042,SCLC-A,SCLC-N_5,SCLC-N_7
RU1083,NSCLC,SCLC-A_3,SCLC-A_1,SCLC-N_4
RU1181,SCLC-A,SCLC-N_6
RU1215,SCLC-A,SCLC-N_2,SCLC-N_3
RU1250,RB1-proficient NSCLC_4,RB1-deficient NSCLC_1,RB1-deficient NSCLC_4,RB1-deficient NSCLC_2
RU1293,NSCLC,SCLC-N_6,SCLC-N_7
RU1303,NSCLC_4,SCLC-A_1,SCLC-A_2
RU1304,SCLC-A,SCLC-N_7,SCLC-N_5
RU1444,NSCLC,SCLC-A_1,SCLC-N_3,UBA52+ SCLC_2
RU1518,NSCLC_1,SCLC-N_8
RU1676,SCLC-A,SCLC-N_7
RU1646,NSCLC,SCLC-N_3
"""
    lines = data.strip().split('\n')
    result_dict = {}
    for line in lines:
        parts = line.split(',')
        key = parts[0]
        values = parts[1:]
        result_dict[key] = values
    samples = ['RU1083', 'RU263', 'RU942', 'RU1444', 'RU1518', 'RU151', 'RU1293', 'RU1303', 'RU581', 'RU831']
    downsample = 1500

    clusters = args.clusters.split("_")

    h5ad_files = [os.path.join(args.path, x + ".no_cc.hvg_2000.090124.h5ad") for x in samples]
    print(f"Creating cellrank probability dataframe with samples: {samples} and downsample: {downsample}")
    df = create_cellrank_probability_df(h5ad_files, 'cell_type_final2', samples, 
                                        result_dict, clusters, downsample=downsample)
    align_obj = CellAlignDTW(df, clusters, 'sample', 'score', 'cell_type')

    print("Aligning data...")
    align_obj.align()

    markers = ["NKX2-1", "NAPSA", "KRT17", "S100A6", "SFTA2", "LGALS3", "SLPI", "AGR2", "TMSB4X", "YAP1",
              "WWTR1", "CHGA", "CHGB", "SYP", "ENO2", "NCAM1", "STMN1", "ASCL1", "NEUROD1", "POU2F3", 
              "INSM1", "MYC", "MYCL", "VIM", "FN1", "RB1", "CDKN2A", "CDKN2B", "CCND1", "CCND2"]
    DE_genes = pd.read_table(args.gene_list)
    DE_genes = DE_genes.iloc[:,0].values
    all_genes = set(markers + DE_genes.tolist())

    print(f"Smoothing expression for genes: {all_genes}")
    summary_df, gene_curves, scores_df = gam_smooth_expression(align_obj.df, 
                                                               all_genes, 
                                                               n_splines=6, lam=3)
    print("Saving results...")
    summary_df.to_csv(args.clusters + '_summary_df.csv')
    gene_curves.to_csv(args.clusters  + "_aggregated_curves.csv")
    scores_df.to_csv(args.clusters  + "_scores_df.csv")

if __name__ == "__main__":
    main()