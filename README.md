# LimmaRP Shiny App #

The app is designed to determine differentially regulated features in data with few replicates and presence of missing values. 
You can test it at [computproteomics.bmb.sdu.dk](http://computproteomics.bmb.sdu.dk)

### Description ###
Detection of differentially regulated features becomes tricky for data with few replicates and high amounts of missing 
values. By combining the commonly used moderated t-test (limma) with rank products based on well-defined null hypotheses, these
features can still be detected with high confidence.

We recommend this method for datasets with a minimum of 1000 features (rows of the table) and at least 3 replicates. The method is based on ratios between features within the replicate, i.e. the tests are carried out on paired tests.

### Literature ###
For details, see Schwämmle, V.; Leon, I. R. and Jensen, O. N. Assessment and improvement of statistical tools for comparative 
proteomics analysis of sparse data sets with few experimental replicates J Proteome Res, 2013, 12, 3874-3883.
[PubMed](http://www.ncbi.nlm.nih.gov/pubmed/23875961)

### Data input ###
Data table (csv file) with optionally row and column names. The values should be log-transformed intensity/abundance values (not ratios) which have been normalized to be comparable. The order of the columns is required to be A1, A2, A3, ..., B1, B2, B3, ..., where 1,2,3 ... are the conditions and A,B,... denote replicates. The tests check for differentially regulated features versus the \"reference\" condition. For each comparison (log-ratio), we plot the histograms of the uncorrected p-values and volcano plots of the false discovery rates (corrected for multiple testing according to Storey JD. A direct approach to false discovery rates. Journal of the Royal Statistical Society. 2002;64:479–498).

### How do I run the program? ###
The easiest way is to install [RStudio](https://www.rstudio.com/) on your computer and run it from there. Alternatively, you can run the files on a [shiny server](https://www.rstudio.com/products/shiny/shiny-server2/) environment.

### Issues ###
For bug reports and general problems, please submit an issue.