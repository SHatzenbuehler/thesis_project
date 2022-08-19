# thesis_project
This repository contains the scripts used for preprocessing and deconvolution analysis for my master thesis project.

The data used here is publicly available from [The Cancer Genome Atlas.](https://portal.gdc.cancer.gov/projects/TCGA-GBM)

The main deconvolution algorithm used here is [MeDeCom (Lutsik et al. 2017)](https://doi.org/10.1186/s13059-017-1182-6) and the repository can be found [here.](https://github.com/lutsik/MeDeCom)


## Non-negtaive Matrix Deconvolution analysis of DNA methylation data for glioblastoma patients (from TCGA).
### Abstract
Glioblastoma multiforme (GBM) is a devastating disease with incidence rates almost as high as mortality and a median survival time of 12 months after diagnosis. Treatment options are very limited and many patients show resistance to the standard approach in chemotherapy. The difficulties faced for new treatment development can be directly tied to tumour heterogeneity between patients as well as intratumoural. Some progress has been made over the years in elucidating this heterogeneity better and there are several subtypes of GBM established in research and clinic. However, lack of characterisation of these subtypes as well as missing information about proportions of cell populations within these subtypes slow the development of subtype-specific treatment approaches.
A way to derive latent cell populations and their proportions from bulk samples is through non-negative matrix factorisation-based deconvolution (NMF). Using DNA methylation samples from The Cancer Genome Atlas, this study investigated the application of NMF in extracting GBM subtype-related information from the data, applying several feature selection methods while doing so.
In this report, distinct latent methylation components are deduced from the data and linked with established GBM subtypes. The proportions of these components are clustered and their usefulness in defining patient groups in combination with subtype classification is investigated. Methodically, this study also highlights differences in variance-based feature selection for NMF compared to other criteria with potential impact on established deconvolution routine applications.


### Contact
For further information on the project or questions:
s.hatzenbuehler@student.maastrichtuniversity.nl
