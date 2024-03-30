# Overview
Repository containing supplemental files and analysis code for 1KGP satellite project.

## scripts:

Folder contains scripts and jupyter notebooks for generating major figures and results from manuscript.

## Supplemental Files:
All Supplemental data can be found in [Zenodo](https://zenodo.org/records/10896290?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImVlYzBiNmQ3LWRjMzQtNGNlNC05ZTBhLTRlOGE0OWFhZjM2YyIsImRhdGEiOnt9LCJyYW5kb20iOiI5OTgwMTYwZjBkMGFiYTI4YjYxZDhiNzVhNWYzYTRkNSJ9.78Ox7KkliHPP02U4iPBnC3hX5Nx8rcR2rZdTe_6X1pQdJQ9asEmlU07zVgAmCFCXpBE0UdbG9u2RhZHIw8MIlQ)

Smaller Supplemental Files can be accessed through this repo, but large files are stored soley in Zenodo.

**Supplemental File 1.** BLAST results of k-mer concatemers against T2T-CHM13-v2.0.

**Supplemental File 2.** Annotations of centromeres, and telomeres of T2T-CHM13-v2.20. Table of abundance of k-mers in annotated regions as numpy file from BLAST hits. Abundance of k-mers across genome in 100kb bins from BLAST hits as .npz files accessible by example:

    import numpy as np
    dense = np.load("filename.npz")
    dense["chr1"]
  
**Supplemental File 3.** Table of pairwise R2 between simple satellites and table of pairwise interspersion OR between simple satellites. Folder containing QQ plots of negative binomial fit of satellite copy number distribution used to qualitatively asses model fit.

**Supplemental File 4.** Materials and results of cenGRM analysis. Boundaries used for centromeric regions of each cenGRM, cenGRMs in GCTA format, and tables with the results of cenGRM GCTA runs. Also provide pdfs of the dendrograms/heatmaps produced from UPGMA clustering of each cenGRM. 

**Supplemental Table 1.** Copy number normalized to 1x depth given GC bias of 126 most abundant satellites analyzed in paper in each individual. Additional columns represent metadata of the individual:

* instrument: sequencer instrument name used to sequence library.
* run: sequencer run of the library.
* flow: flowcell ID of the ibrary.
* pop: 1,000 Genomes Project population ID.
* superpop: 1,000 Genomes Project superpopulation ID.
* reads: average autosomal read depth of the library.

**Supplemental Table 2.** Copy number normalized to 1x depth given GC bias of the top 126 most abundant satellites analyzed in paper in each individual of the 1KGP, plus estimates of the same satellites in CHM13 short-read libraries subsampled from 18x-0.5x, 18x depth simulated library of the T2T-CHM13v2.0 assembly analyzed using k-Seek, and Tandem Repeat Finder results of the T2T-CHM13v2.0 asembly [(Hoyt 2022)](DOI: 10.1126/science.abk3112).

