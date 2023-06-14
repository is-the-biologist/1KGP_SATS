Scripts that are used to perform analysis of k-seek data. Raw data output from k-seek is too large for github, but easily reproducible by following k-seek manual and applying to 1KGP fastqs. The jupyter notebooks included in "git_jupyter" are used to generate the main figures and results of this manuscript. Please use conda (https://docs.conda.io/en/latest/) to run the notebooks and activate the environment.

# Activate environments:

To run the Jupyter notebooks please use the environment ".yml" files in Conda first. This will ensure the notebooks have the same dependencies as on publication. To do this within the "git_jupyter" directory and then:

    conda env create -f HumanGenomeSatVar_gitcode.1.yml
    conda env create -f HumanGenomeSatVar_gitcode.2.yml

Then activate the environments before runing notebook HumanGenomeSatVar_gitcode.1.ibynb and HumanGenomeSatVar_gitcode.2.ibynb, respectively:
    
    conda activate HumanGenomeSatVar_gitcode.1

or
 
    conda activate HumanGenomeSatVar_gitcode.1    
 

to activate 

# Jupyter Notebooks:

# Python scripts

* interspersion.py: Script contains algorithm for counting interspersed satellites based on the read-pair output from k-seek.
* normKseek.py: Describes the algorithm used for normalizing raw k-mer counts data by GC and read depth using a GC bin table.
