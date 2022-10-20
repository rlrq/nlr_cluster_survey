# nlr_cluster_survey

~~(Note that the files currently in the repo reproduce the analysis in the manuscript before major revisions. Final scripts are to be uploaded in due time.) (Now with CNVnator analysis for A. lyrata.)~~ The scripts should be final, but human error may mean that something could have been mistakenly left out (or in) during the reorganisation of code for readability.

pipeline_cluster_survey.sh details the analysis pipeline up to but not including visualisation steps. Miscellaneous scripts required are included in the 'scripts' folder, and the directory structure should be preserved.

Paths to update before execution
- Variables at the top of pipeline_cluster_survey.sh
- Variables at the top of plot_cluster_survey.R
- Paths to FASTA files in scripts/get_seqs/scripts/get_seqs_functions.py

Manual steps
- The NB-ARC tree should be manually rooted by assigning the division between TIR and non-TIR NLRs as the most basal before proceeding with the "assign genes" step.


The pipeline was executed with
- Ubuntu 18.04.3
  - bash 4.4.20
- Python 3.6.9
- R 3.6.1

Third party command line tools
- bedtools v2.26.0
- BLAST+ 2.6.0
- mafft v7.427
- FastTree 2.1
- CNVnator
  - Requires: Yeppp!

Third party Python libraries
- Biopython 1.73

Third party R packages
- bedr v1.0.7
- tidyverse 1.2.1
- data.table 1.12.8

Additional third party R packages for plotting and visualisation (script not yet uploaded)
- grid 3.6.2
- gridExtra 2.3
- ggpubr 0.2.4
- wesanderson 0.3.6
- RColorBrewer 1.1.2
- Biostrings 2.52.0
- rlist 0.4.6.1
- ape 5.3
- ggsignif 0.6.0
- ggpmisc 0.3.3
- geiger 2.0.6.3
- scales 1.1.0
- ggtree 1.16.6
- ggridges 0.5.2
