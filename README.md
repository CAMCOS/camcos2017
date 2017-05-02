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
