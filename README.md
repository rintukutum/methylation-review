# methylation-review

#### Steps:

```{shell}
# terminal
# Find combination
Rscript A-find-combination-methylation.R

# Find variability
Rscript B-variability-per-gene.R
```

#### Notes
Please contact the authors for `raw-filtered-data.csv` to reproduce the results. Please install the pre-requisites before running the scripts.
```{R}
install.packages(
  c(
    "randomForest",
    "Boruta",
    "plyr",
    "reshape",
    "ggplot2"
  )
)
```
