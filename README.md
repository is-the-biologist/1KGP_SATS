# 1KGP_SATS
Repository containing supplemental files and analysis code for 1KGP satellite project.


**Supplemental Files:**

**Supplemental File 1.** BLAST results of k-mer concatemers against T2T-CHM13-v2.0. (Files are stored in google drive link at end of section due to Github size limit).

**Supplemental File 2.** Annotations of centromeres, and telomeres of T2T-CHM13-v2.20. Table of abundance of k-mers in annotated regions as numpy file from BLAST hits. Abundance of k-mers across genome in 100kb bins from BLAST hits as .npz files accessible by example:

    import numpy as np
    dense = np.load("filename.npz")
    dense["chr1"]
  
**Supplemental File 3.** Table of pairwise R2 between simple satellites and table of pairwise interspersion OR between simple satellites. Folder containing QQ plots of negative binomial fit of satellite copy number distribution used to qualitatively asses model fit.

**Supplemental File 4.** Materials and results of cenGRM analysis. Boundaries used for centromeric regions of each cenGRM, cenGRMs in GCTA format, and tables with the results of cenGRM GCTA runs. Also provide pdfs of the dendrograms/heatmaps produced from UPGMA clustering of each cenGRM. (Files are stored in google drive link at end of section due to Github size limit).

**Supplemental Table 1.** Copy number normalized to 1x depth given GC bias of 126 most satellites analyzed in paper in each individual. Additional columns represent metadata of the individual:

* instrument: sequencer instrument name used to sequence library.
* run: sequencer run of the library.
* flow: flowcell ID of the ibrary.
* pop: 1,000 Genomes Project population ID.
* superpop: 1,000 Genomes Project superpopulation ID.
* reads: average autosomal read depth of the library.

**Google Drive link for Supplemental File 1 and 4.** https://drive.google.com/drive/folders/1g5dvGMHUw4lazEPD48clsNu5STNPPR3G?usp=sharing
