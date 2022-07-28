--- 
title: "DNA methyltransferases 3A and 3B target specific sequences during mouse gastrulation"
author: "Aviezer Lifshitz"
date: "2022-07-28"
knit: "bookdown::render_book"
site: bookdown::gitbook_site
description: "Analysis and code for the paper: DNA methyltransferases 3A and 3B target specific sequences during mouse gastrulation"
url: 'tanaylab.github.io/Dnmt3ab_MEEB'
github-repo: tanaylab/Dnmt3ab_MEEB
---

# DNA methyltransferases 3A and 3B target specific sequences during mouse gastrulation

Following is the code that generates the figures for the methylation part of Mukamel et al. paper on DNA methyltransferases 3A and 3B activity during gastrulation. The code is splitted to jupyter notebooks that can be found at: https://github.com/tanaylab/Dnmt3ab_MEEB

The RNA part can be found at that repository under the `rna` directory.

## Run the notebooks

Prior to any analysis, after cloning the repository, please download first the necessary data by running (in the root directory of the cloned repository):


```bash
R -e "source('scripts/download_data.R'); download_full_data()"
```

## Download the UMI matrices 

The UMI matrices can be downloaded by running: 

```bash
R -e "source('scripts/download_data.R'); download_raw_umi_tables()"
```

## Notebook order

- Methylation-trends.ipynb
- Sequence-model.ipynb
- Sequence-model-validation.ipynb
- Strand-specific-models.ipynb
- Yagi_et_al.ipynb
- Single-cell-methylation.ipynb
- Day6-differential-methylation.ipynb
- QC-bulk-meth.ipynb
- QC-sc-meth.ipynb
