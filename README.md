# 1KGP_SATS
Repository of satellite variation in 1KGP data.


Supplemental Files:

Supplemental File 1. BLAST results of k-mer concatemers against T2T-CHM13-v2.0.

Supplemental File 2. Annotations of centromeres, and telomeres of T2T-CHM13-v2.20. Table of abundance of k-mers in annotated regions as numpy file from BLAST hits. Abundance of k-mers across genome in 100kb bins from BLAST hits as .npz files accessible by example:

  import numpy as np
  dense = np.load("filename.npz")
  dense["chr1"]
  
Supplemental File 3. Table of pairwise R2 between simple satellites and table of pairwise interspersion OR between simple satellites. Folder containing QQ plots of negative binomial fit of satellite copy number distribution used to qualitatively asses model fit.

Supplemental File 4. Materials and results of cenGRM analysis. Boundaries used for centromeric regions of each cenGRM, cenGRMs in GCTA format, and tables with the results of cenGRM GCTA runs. Also provide pdfs of the dendrograms/heatmaps produced from UPGMA clustering of each cenGRM.

