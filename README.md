```{bash}
# Download scripts have to be placed in 'download_scipts' directory.
# Callable from the dataset's root folder, e.g like so:

# use the download+normalize+remove batch effect script on the samples
# defined in samples/whole_blood.R
Rscript download_scripts/dl_normalize_rmbatch.R samples/whole_blood.R
```
