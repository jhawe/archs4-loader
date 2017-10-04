## Overview
This small collection of scripts makes use of the [ARCHS<sup>4</sup> resource](http://amp.pharm.mssm.edu/archs4/)<sup>[1](https://www.biorxiv.org/content/early/2017/09/15/189092)</sup> ([What is ARCHS<sup>4</sup>?](https://github.com/MaayanLab/archs4-extension#what-is-archs4)).
The scripts allows for an easy (rather naive) definition of samples which are to be extracted. Additionally,
a definition of samples can be obtained from the [ARCHS<sup>4</sup> data portal](http://amp.pharm.mssm.edu/archs4/data.html).
The download script is adjusted from the scripts which are provided by the ARCHS<sup>4</sup> website, allowing for optional normalization and batch effect removal.
Furthermore, simple diagnostic plots are created for the extracted data.

## Defining samples
You can create definitions of which samples to be used for extracting expression data (i.e., files
located under `./sample_definitions/`) by using the `./scripts/getSamples.R`.
The given keywords are currently searched independently (in an 'either or' fashion) in the 'tissue' meta data field only.
If you want more keywords connected by 'AND', enclose them with quotes like in the example below.

### Example
```{bash}
Rscript scripts/getSamples.R "Whole blood"
Rscript scripts/getSamples.R heart
```

## Extracting data
The download script is based on the human_matrix_download.h5 file located in the root directory.
The *.h5 file is an archive file which contains all data and meta data relevant to get the 
gene expression data from the ARCHS<sup>4</sup> repository. If the file does not yet exist in the repositories root directors (it is NOT provided in the repository at the moment) then it gets downloade from the ARCHS<sup>4</sup> website once. 
The download script (`./scripts/download.R`) checks whether the file already exists and then uses the file to extract the expression data for the specified samples. The script extracts the (normalized) expression data for the specified samples
and creates some basic diagnostic plots in a subfolded named according to the sample-definition script.

### Available arguments

|Flag|Meaning|
|-----|------|
|`--samples`|The sample definition file which can be sourced in R. Compare files in 'sample_definitions/'|
|`--normalize`|Whether to perform log-transform and quantile normalization|
|`--sva`|Whether to use SVA based batch effect removal|
|`--peer`|Whether to use PEER based batch effect removal|

## Example
```{bash}
# download whole_blood samples
Rscript scripts/download.R --normalize --sva --samples=sample_definitions/whole_blood.R
# download whole_blood samples (raw counts)
Rscript scripts/download.R --samples=sample_definitions/whole_blood.R
# download heart samples using peer
Rscript scripts/download.R --normalize --peer --samples=sample_definitions/heart.R
```

## TODOs

* versioning of h5 file (e.g. by month/year)
* check for ChIP-seq data as well
* check batch effect removal of normalized expression data
