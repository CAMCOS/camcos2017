# CAMCOS Spring 2017

### Install

Install the release version of `devtools` from CRAN with `install.packages("devtools")`. 
Then run the following line of code to install the `camcos2017`

```R
devtools::install_github("CAMCOS/camcos2017")
```

### Dependencies

* Data manipulation:
    * `dplyr`
* Clustering
    * `RSpectra`
    * `Matirx`
* Vizualization
    * `Rtsne`
    * `ggplot2`

### Load Data

The full 20newsgroups dataset (test and train) is stored in the data object: `newsgroups`. To load the data:

```R
data(newsgroups)
```

To select a subset of the full dataset, binarize, or append the vocabulary:

```R
## Filter by subject area
sub0 <- news_subset(newsgroups, filter = c("comp","rec"), binary = TRUE, vocab = TRUE)

## Filter by label numbers
sub1 <- news_subset(newsgroups, filter = 2:6, binary = TRUE, vocab = TRUE)
```

### Perform Clustering

```R
library(magrittr)
x <- colweights(data = sub0[,1:3], binary = FALSE, weightfunction = "IDF") %>%
  similarity("correlation") %>%
  clustering("DiffusionMap", k = 6, t = 0.5) %>%
  clustercheck(labels.6, k = 6)
cat("Accuracy: ", round(x[[3]], 2), "%", sep = "", fill = TRUE)
```
