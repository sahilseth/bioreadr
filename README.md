ngsflows
------------------------------

A collection of NGS related functions which wrap several well known tools.

One may string together a pipeline using these functions.


## Downloading resource bundles.

Reference libraries are available from Illumina and [Broad Insititute](https://www.broadinstitute.org/gatk/guide/article.php?id=1215).

We have a in built function to get genomes from Illumina (example below), but for Human genetics its better to use the bundle from Broad Institute.

One may run the following commands to get the bundle from Illumina

```
library(ngsflows)
fetch_genomes(species="Homo_sapiens")

################################################
Available Sources:
################################################

Ensembl
NCBI
UCSC

Example:
fetch_genomes species=Homo_sapiens src=NCBI build=build37.2

```


```
wget ... and download the reference data and files
## wget download the 
## test R
Rscript --version
## get runflow from github, name of the flow, samplesheet, parameterfile
## wget runflow
```


```r
install.packages('devtools')
require(devtools)
install_github(repo = 'ngsflows', username = 'sahilseth')
```
