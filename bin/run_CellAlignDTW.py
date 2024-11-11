import os
import argparse
import pandas as pd
from CellAlignDTW import CellAlignDTW, gam_smooth_expression
from CellAlignDTW.pp import stratified_downsample

def main():
    parser = argparse.ArgumentParser(description="Run CellAlignDTW analysis")
    parser.add_argument("--path", required=True, help="Path to the directory containing the files")
    parser.add_argument("--outpath", required=True, help="Path to the output directory")
    parser.add_argument("--clusters", required=True, help="Comma-separated list of clusters to analyze")
    #parser.add_argument("--patients", required=True, help="Comma-separated list of patients to analyze")
    parser.add_argument("--gene_list", required=True, help="Path to the file containing the list of DE genes")
    parser.add_argument("--n_splines", required=True, help="n splines for GAM", type=int)
    parser.add_argument("--lam", required=True, help="lambda for GAM", type=float)
    args = parser.parse_args()

    print(f"Starting CellAlignDTW analysis with path: {args.path}, clusters: {args.clusters}, and gene list: {args.gene_list}")

    downsample = 1500

    clusters = args.clusters.split("_")
    
    df = pd.read_csv(os.path.join(args.path, f'cellrank.{args.clusters}.csv'))
    df = df.groupby('sample').apply(lambda group: stratified_downsample(group, 'score', downsample)).reset_index(drop=True).drop(columns=['Unnamed: 0'])
    align_obj = CellAlignDTW(df, clusters, 'sample', 'score', 'cell_id', 'cell_type')

    print("Aligning data...")
    align_obj.align()

    all_genes = ["NKX2-1", "NAPSA", "KRT17", "S100A6", "SFTA2", "LGALS3", "SLPI", "AGR2", "TMSB4X", "YAP1",
              "WWTR1", "CHGA", "CHGB", "SYP", "ENO2", "NCAM1", "STMN1", "ASCL1", "NEUROD1", "POU2F3", 
              "INSM1", "MYC", "MYCL", "VIM", "FN1", "RB1", "CDKN2A", "CDKN2B", "CCND1", "CCND2", "PHOX2B"]

    DE_genes = pd.read_table(args.gene_list, sep=None, engine='python')
    if DE_genes.shape[1] > 1:
        DE_genes = DE_genes.iloc[:,0].values
        all_genes = set(all_genes + DE_genes.tolist())
    
    print(f"Smoothing expression for genes: {all_genes}")
    summary_df, gene_curves, scores_df = gam_smooth_expression(align_obj.df, 
                                                               all_genes, 
                                                               n_splines=args.n_splines, 
                                                               lam=args.lam)
    print("Saving results...")
    summary_df.to_csv(os.path.join(args.outpath, args.clusters + '_summary_df.csv'))
    gene_curves.to_csv(os.path.join(args.outpath, args.clusters  + "_aggregated_curves.csv"))
    scores_df.to_csv(os.path.join(args.outpath, args.clusters  + "_scores_df.csv"))

if __name__ == "__main__":
    main()