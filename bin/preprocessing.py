import os
import argparse
from scipy.io import mmwrite
import pandas as pd
from CellAlignDTW.pp import create_cellrank_probability_df

def main():
    parser = argparse.ArgumentParser(description="Preprocessing for CellAlignDTW analysis")
    parser.add_argument("--path", required=True, help="Path to the directory containing the h5ad files")
    parser.add_argument("--clusters", required=True, help="Clusters to analyze (concatenated by underscores)")
    parser.add_argument("--outpath", required=True, help="Path to the output directory")
    args = parser.parse_args()

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
    clusters = args.clusters.split("_")
    h5ad_files = [os.path.join(args.path, x + ".no_cc.hvg_2000.090124.h5ad") for x in samples]
    combined_adata, df = create_cellrank_probability_df(h5ad_files, 'cell_type_final2', samples, 
                                        result_dict, clusters)
    print(f"Pseudobulking for cell types: {clusters}")
    # Create a DataFrame from the counts layer
    counts_df = pd.DataFrame(combined_adata.layers['counts'], 
                                                index=combined_adata.obs_names, 
                                                columns=combined_adata.var_names)
    counts_df['sample'] = df['sample'].values
    counts_df['cell_type'] = df['cell_type'].values
    summed_counts = counts_df.groupby(['sample', 'cell_type']).sum().reset_index()
    print(summed_counts.shape)

    print("Saving results...")
    summed_counts.to_csv(os.path.join(args.outpath, f'counts.{args.clusters}.csv'))
    df.to_csv(os.path.join(args.outpath, f'cellrank.{args.clusters}.csv'))

if __name__ == "__main__":
    main()