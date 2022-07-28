# DNA methyltransferases 3A and 3B target specific sequences during mouse gastrulation

Files in this directory are supporting Mukamel et al. work on Dnmt3a/b activity during gastrulation.

Scripts in the `rna` directory can reproduce metacell derivation and analysis generating figure panels 1-3. The metacell models used in the papers are also provided.

DNA methylation data is provided in a processed form to ease up analysis, and the scripts that reproduce it are provided as jupyter notebooks at `methylation/analysis` directory. 

The processed files can be downloaded by running: 

```bash
R -e "source('scripts/download_data.R'); download_full_data()"
```

Or if you only want the UMI matrices: 

```bash
R -e "source('scripts/download_data.R'); download_raw_umi_tables()"
```

For both scRNA-seq and 5mC data, raw data has been deposited at GEO (GSE199806).

A running [MCView](https://tanaylab.github.io/MCView/) app is live at: 

https://tanaylab.weizmann.ac.il/Dnmt3ab_MEEB

## Dependencies

R packages (available from CRAN): 

- tidyverse
- glue
- here
- patchwork
- ggseqlogo
- tgstat
- tglkmeans
- xgboost
- caret
- glmnet
- devtools
- ggrepel
- qvalue
- scales
- ggsci
- writexl
- ggforce
- matrixStats
- ggtext
- data.table
- ggridges
- princurve
- zoo
- skimr
- broom
- viridis
- tidyr
- ggplot2 
- dplyr 
- forcats
- purrr
- plyr
- scattermore
- readr
- ggpubr
- cowplot
- rlang
- doMC
- grDevices
- pheatmap

R packages (available from github):

- tanaylab/metacell 
- tanaylab/tgutil 
- tanaylab/tgppt 
- tanaylab/misha 
- tanaylab/tgconfig 
- tanaylab/misha.ext
- tanaylab/gpatterns
- tanaylab/sc5mc # note that this package requires intel MKL libraries to be installed
- tanaylab/tgcg

BioConductor R package:

- ComplexHeatmap