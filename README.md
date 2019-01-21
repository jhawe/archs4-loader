## Overview
This small collection of scripts makes use of the [ARCHS<sup>4</sup> resource](http://amp.pharm.mssm.edu/archs4/)<sup>[1](https://www.biorxiv.org/content/early/2017/09/15/189092)</sup> ([What is ARCHS<sup>4</sup>?](https://github.com/MaayanLab/archs4-extension#what-is-archs4)).
The scripts allow for an easy (but rather naive) definition of samples which are to be extracted. Additionally,
a definition of samples can be obtained from the [ARCHS<sup>4</sup> data portal](http://amp.pharm.mssm.edu/archs4/data.html).
The download script is adjusted from the scripts which are provided by the ARCHS<sup>4</sup> website, allowing for optional normalization and batch effect removal.
Furthermore, simple diagnostic plots are created for the extracted data.

## Creating a sample sheet
To create a new sample/design sheet, simply call the script 'create_sample_sheet.R' with no arguments.
This will 

- download the main h5 file
- extract information from ALL samples contained in this file

```{bash}
Rscript scripts/create_sample_sheet.R
```

## Defining samples
You can create definitions of which samples to be used for extracting expression data (i.e., files
located under `./sample_definitions/`) by using the `./scripts/getSamples.R` R script.
The given keywords are currently searched independently (in an 'either or' fashion) in the 'tissue' meta data field only.
If you want more keywords connected by 'AND', enclose them with quotes like in the example below.

### Example
```{bash}
Rscript scripts/getSamples.R "Whole blood"
Rscript scripts/getSamples.R heart
Rscript scripts/getSamples.R "Skeletal Muscle"
```

## Extracting data
The download script is based on the human_matrix_download.h5 file located in (and independently downloaded to) the root directory of the project.
The \*.h5 file is an archive file which contains all data and meta data relevant to get the 
gene expression data from the ARCHS<sup>4</sup> repository. If the file does not yet exist in the repositories root directors (it is NOT provided in the repository at the moment) then it gets downloade from the ARCHS<sup>4</sup> website once. 
The download script (`./scripts/download.R`) (and the sample sheet creation script) check whether the file already exists and then use the file to extract the relevant sample design or expression data. The download script extracts the (normalized) expression data for the specified samples and creates some basic diagnostic plots in a subfolded named according to the sample-definition script.

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
