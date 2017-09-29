The scripts are based on the human_matrix_download.h5 file based in the root directory.
This is an archive file which contains all data and meta data relevant to get the 
gene expression data. The download scripts should check whether the file already exists 
(otherwise download it) and then use the file to extract the expression data. 
Download scripts have to be placed in 'download_scripts' directory and be callable from 
the dataset's root folder, e.g like indicated below.

```{bash}
# use the download+normalize+remove batch effect script on the samples
# defined in samples/whole_blood.R
Rscript download_scripts/dl_normalize_rmbatch.R samples/whole_blood.R
```

TODOs: 

* versioning of h5 file (e.g. by month/year)
* parse metadata for additional experiment information and create "design" file for dls
