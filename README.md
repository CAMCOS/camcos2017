# CAMCOS Spring 2017

### Install

Install the release version of `devtools` from CRAN with `install.packages("devtools")`. 
Then run the following line of code to install the `camcos2017`

```R
devtools::install_github("CAMCOS/camcos2017")
```

### Dependencies

The function `rtsne_viz` requires the package `Rtsne`, which should be automatically installed when the `camcos2017` package is installed. 

The function `similarity` requires the package `proxy`. This is **not** automatically installed with this package because it requires version 3.3.2 or higher of R. If you have the appropriate version of R, install the `proxy` package manually with `install.packages("proxy")` before using the function.

### Load Data

The data subset containing one cluster from each of the main cluster groups is included in the package. 
The first line of code below will load two objects into the R environment: `cluster6_counts` and `cluster6_labels`. 
The next line binarizes the counts by setting all values to 1.

```R
data(cluster6)
cluster6_counts$value <- 1
```

Note: after extracting the six clusters from the original training data set, only two processing steps were performed:

1. Remove any columns that have a column sum of 0
2. Remap the original cluster labels to 1-6: {2: 1, 7:2, 8:3, 12:4, 19:5, 20:6}
