# nlr_cluster_survey

Paths to update before execution
- Variables at the top of pipeline_cluster_survey.sh
- FASTA files in scripts/get_seqs/scripts/get_seqs_functions.py

Manual steps
- The NB-ARC tree should be manually rooted by assigning the division between TIR and non-TIR NLRs as the most basal before proceeding with the "assign genes" step.


The pipeline was executed with
- Ubuntu 18.04.3
  - bash 4.4.20
- Python 3.6.9
- R 3.6.1

Third party command line tools
- bedtools v2.26.0
- blast 2.6.0
- mafft v7.427
- FastTree 2.1

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
