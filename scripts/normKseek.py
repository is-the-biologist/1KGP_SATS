import pandas as pd
import numpy as np
import os
import sys

def reformatF(file):
    k_df = pd.read_csv(file, index_col=0)
    k_df.columns = [k.replace(".", "/") for k in k_df.columns]
    k_df['sample'] = k_df.index
    k_df.to_csv(file.replace('.csv', '.compiled'), index=None, sep="\t")

def normalize(fname, outfname, set_bounds=False):

    tandems = pd.read_csv(fname, sep="\t")
    gc_bins = np.asarray([ ((b+1)/100) for b in range(100)])
    if set_bounds: #bounds are not used in any run -- ignor this param
        get_gc_bins = lambda gc: min(gc_bins[  min(max(float(gc), 0.015), 0.79) <= gc_bins]) #gets gc bin for a given % GC and also bounds GC estimate to reasonable bins.
        #bounding is optional and will be parameter.
    else:
        get_gc_bins = lambda gc: min(gc_bins[float(gc) <= gc_bins])
    #get GC content of each monomer and implement a function to trim the edges of the GC estimate bins bc I'm not confident of GC estimation at those edges:
    GC_monomers = [  get_gc_bins( (np.sum(np.asarray(list(k.split("/")[0])) == "G") + np.sum(np.asarray(list(k.split("/")[0])) == "C") )/ len(k.split("/")[0]) ) for k in tandems.loc[:, ~tandems.columns.isin(['total_bp', 'sample'])].columns if np.sum(np.asarray(list(k)) == "N") == 0]

    cols = [k for k in tandems.loc[:, ~tandems.columns.isin(['total_bp', 'sample'])].columns if np.sum(np.asarray(list(k)) == "N") == 0]

    correction_factor_matrix = np.full(fill_value=np.nan, shape=(tandems.shape[0], len(cols)))
    ind = 0
    sample_order = tandems['sample'].values
    for f in sample_order:
        #append all correction factors to matrix
        depth_summary = pd.read_csv( os.path.join("/fs/cbsuclarkfs1/storage/30x_1K_genomes/gc/GC_EST_UNIQ_2022_05_18", f+".gc.depth.tsv"), sep="\t" , header=None)

        #calculate depth | %GC using bins
        med_value = np.median(depth_summary[3].values)
        #adding conditional argument such that any bins that had < 10000 reads compose them are to be discarded and instead use the median value
        #this is because the edges of the estimates become spikey and unstable as the # of reads that make up a bin
        #is quite small
        gc_depth = {d[0]: (d[3] if d[1] >= 10000 else med_value) for d in depth_summary.values} #assign bins to dict
        #calculate the normalization factor based on GC content
        norm_factor = [gc_depth[GC] for GC in GC_monomers]
        correction_factor_matrix[ind,:] = norm_factor
        ind += 1

    #match order of sample sheet in matrix
    copy_number = tandems.set_index('sample').loc[sample_order][cols].values #probably depreceated now
    normalized_copies = copy_number / correction_factor_matrix

    norm_tandem_df = pd.DataFrame(data=normalized_copies, columns=cols)

   #set index of new dataframe to sample_order which is the correct labels
    norm_tandem_df["sample"] = sample_order

    norm_tandem_df.to_csv(outfname, index=False)


if __name__ == '__main__':
    arg = sys.argv[1]

    fname_comp = arg.replace('.csv', '.compiled')
    fname_gc = fname_comp.replace("RAW", "GC_CORRECTED")


    reformatF(file=f"/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/RUN_2022_06_28/final_tables/{arg}")
    normalize(fname=f"/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/RUN_2022_06_28/final_tables/{arg.replace('.csv', '.compiled')}",
              outfname=f"/fs/cbsuclarkfs1/storage/is372/human_kmer/full_kseek_run/RUN_2022_06_28/final_tables/{fname_gc}")