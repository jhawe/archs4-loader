# Overview
This small collection of scripts makes use of the [ARCHS<sup>4</sup> resource](http://amp.pharm.mssm.edu/archs4/)<sup>[1](https://www.biorxiv.org/content/early/2017/09/15/189092)</sup> ([What is ARCHS<sup>4</sup>?](https://github.com/MaayanLab/archs4-extension#what-is-archs4)).
The scripts allow for an easy (but rather naive) definition of samples which are to be extracted. Additionally,
a definition of samples can be obtained from the [ARCHS<sup>4</sup> data portal](http://amp.pharm.mssm.edu/archs4/data.html).
The download script is adjusted from the scripts which are provided by the ARCHS<sup>4</sup> website, allowing for optional normalization and batch effect removal.
Furthermore, simple diagnostic plots are created for the extracted data.

# Snakemake workflow 
## Data download and processing
We use a [Snakemake](https://snakemake.readthedocs.io/) based approach to have access to a fully automated workflow.
To only download the data, as well as perform batch effect removal, the pipeline can be called by executing the following line of code in the root directory of the project:

```{bash}
snakemake results/downloads/{your_keywords}/design.tsv
```

This will obtain all samples from ARCHS<sup>4</sup> which have the specified keywords (separated by "_")
annotated in their tissue meta-data field.
As of now, the downloaded expression data will be automatically normalized using [ComBat](https://www.bu.edu/jlab/wp-assets/ComBat/Abstract.html)

## Data exploration
After downloading and processing the data, we now can go one step further and create a basic
overview. This is implemented as an Rmarkdown and is now integrated in the snakemake workflow:

```{bash}
snakemake results/downloads/{your_keywords}/summary.html
```
## Note
Snakemake version 5.2.2 or greater is required for the Rmarkdown to render properly.
This version contains a bugfix crucial for successfully loading the data

# TODOs
- allow parametrization (e.g. type of normalization, raw data saving etc)
