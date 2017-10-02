## Overview
This small collection of scripts makes use of the [ARCHS4 resource](http://amp.pharm.mssm.edu/archs4/).
The scripts allows for an easy (rather naive) definition of samples which are to be extracted. Additionally
a definition of samples can be obtained from the [ARCHS4 data portal](http://amp.pharm.mssm.edu/archs4/data.html).
The download script is adjusted from the scripts which are provided by the ARCHS4 website, allowing for optional normalization and batch effect removal.
Additionally, simple diagnostic plots are created for the extracted data.

## Defining samples
You can create definitions of which samples to be used for extracting expression data (i.e., files
located under ./sample_definitions/) by using the "getSamples.R" script in the scripts/ subfolder.
The given keywords are currently searched independently (in an 'either or' fashion) in the 'tissue' meta data field only.
If you want more keywords connected by 'AND', enclose them with quotes like in the example below.

```{bash}
Rscript scripts/getSamples.R "Whole blood"
```

## Extracting data
The download scripts are based on the human_matrix_download.h5 file located in the root directory.
This is an archive file which contains all data and meta data relevant to get the 
gene expression data from the archs4 repository. If the file does not yet exist in the repositories root directors (it is NOT provided in the repository at the moment) then it gets downloade from the ARCHS4 website once. 
The download script checks whether the file already exists and then uses the file to extract the expression data for 
the specified samples. The script extracts the (normalized) expression data for the specified samples
and creates some basic diagnostic plots in a subfolded named according to the sample-definition script.

```{bash}
# use the download+normalize+remove batch effect script on the samples
# defined in samples/whole_blood.R
Rscript scripts/download.R --normalize --rmbatch --samples=sample_definitions/whole_blood.R
```

## TODOs

* versioning of h5 file (e.g. by month/year)
* check for ChIP-seq data as well
* check batch effect removal of normalized expression data
