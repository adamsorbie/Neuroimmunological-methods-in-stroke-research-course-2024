# Gut_brain_axis_course-Bioinformatic_analysis



Data analysis part of the 2023 Gut_brain_axis_course-Bioinformatic_analysis - Thursday 29/06/23

## Intro-to-R 

Since some of you may not have any experience using R, we will very quickly go through some key concepts to help you understand better what we will doing in the analysis section later. 

### Mathematical operations 

Type and test these operations in the console in Rstudio 

Addition
```1 + 1```

Subtraction 
```5 - 4``` 

Multiplication 
```2 * 10```

Division
```9 / 3```

Modulo 
```20 %/% 6```

Square
```4^2 or 4**2``` 

### Comparison

Greater than
```>``` 

Less than 
```<``` 

Equals 
```==```

Does not equal 
```!=```

### Basic types 

There are 6 basic data types in R: logical, numeric, integer, character, raw and complex. Today (and in general) only the first four are important. 

```logical``` Can only have two values ```TRUE``` or ```FALSE``` (Shorthand:```T```,```F```). 

```numeric``` All numbers with or without decimal ```100.34```

```integer``` Whole numbers ```123L``` A number appended with L suffix denotes an integer in R

```character``` Character or string values. Enclosed in either single 'microbe' or double "stroke" quotes. 

### Data structures 

Four most commonly used data structures 

#### 1D
Vectors - ordered collection of data of a given length and all of the same type:
``` r
vector1 <- c(1,3,7,8,9)
vector2 <- c("A", "B", "C", "D")
print(vector2)
```

Lists - ordered collection of data/objects of a given length. Can be composed of heterogenous types. 
``` r
list1 <- list(vector1, vector2)
list2 <- list(1,2,3)
print(list1)
```

#### 2D 

Matrices - 2D data structure (rows x columns) of the same type
``` r
m <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3)
print(m)
```

Dataframes - 2D heterogenous data structure (rows x columns)
``` r
df <- data.frame(bacteria = c("E. coli", "L. reuteri", "C. difficile"),
                 sample1 = c(10, 27, 61),
                 sample2 = c(9, 200, 43) )
print(df) 
```

### Variable assignment 

Variables are assigned in R using the ```<-``` operator. A variable can be though of a a "box" you store a value in for later use.

They can be used to store a value which you want to be constant throughout your code e.g. a file name 
``` r
x <- "Bacteria_data.txt" 
``` 
Or to save time by storing long strings/numbers 

``` r
y <- 127575738398292929
```

Objects such as vectors and dataframes can also be stored in variables as you saw above. ```df``` and ```vector2``` are both variables. 

### Control flow

As well as writing individual statements like we have done above we can also use logic to control the flow of our code, allowing us to repeat bits of code or only run if a given condition is met. 

for loops allow us to repeat code a specified number of times e.g.

``` r
for (i in vector2){
  print(i)
 }
```

This prints out each element in the ```vector2```variable we defined earlier. Not particularly interesting though.. 
We will discover a more interesting use case later. 

if/else statements allow us to control the flow of our code better: 

``` r
x <- 100
y <- 10

if (x > y) {
  print("x is higher")
} else if (x < y) {
  print("y is higher")
}
```
Change the values of x and y and see how the output changes. 

### Functions 

Functions are used to abstract repetitive code in something which is reusable and modular, making code easier to read and debug. 

Using the above if else example we could create a function called ```greater_than``` which tells us which value is highest in a more modular way

``` r
greater_than <- function(x, y) {
  if (x > y) {
    print("x is higher")
  } else if (x < y) {
    print("y is higher")
  }
}
```

``` r
greater_than(x=100, y=10)
greater_than(x=1, y=17)
```

More interesting functions might perform a calculation for us and return the value 

``` r
normalise <- function(x) {
  x_norm <- t(100 * t(x) / colSums(x))
  return(x_norm)
}
```
Let's run this function on the matrix we created earlier to see what it does. 

``` r
normalise(m)
```

So far we have covered built-in functions and custom functions, but R has a huge 
open source library of packages called CRAN, as well as a specific repository
for bioinformatics, called bioconductor. 

### Loading external function 

When working with R, at some point you may find yourself requiring a specific 
function to perform some kind of analysis and plot data etc. It's usually the 
case that someone has already written a function for whatever you want to do

You can install packages to access functions written by others. 
External packages can be downloaded using the ```install.packages``` function 
and loaded using the ```library``` function. Somewhat annoyingly, bioconductor
packages are installed a little differently, and require you to first install
the "BiocManager" package and then install using ```biocManager::install```. They
can however be loaded as normal. 

If you need help with a function you can also type ?functionname in the console e.g. ?log10 and the help for that function will show up, detailing what the function does, what inputs it expects and what value(s) it returns. 

Once you get into writing your own functions, it's good practice to store them
in a separate R file and import them using the ```source``` function. 

### Data wrangling 

A few key concepts on loading and manipulating data. 

Reading data 
``` r
# load tidyverse
library(tidyverse)
# read in the palmer penguins csv file 
penguins <- read_csv("data/penguins.csv")
```

Manipulating data 

Selecting columns 

Base R 

``` r
bill_length <- penguins$bill_length_mm
print(bill_length)
# alternative
bill_length <- penguins["bill_length_mm"]
print(bill_length)
```
Tidyverse
``` r
# select certain columns 
bill_length <- select(penguins, bill_length_mm) 
print(bill_length)
# using pipe operator
bill_length <- penguins %>% 
    select(bill_length_mm)
print(bill_length)
```

Filtering rows 

Base R 
``` r
gentoo  <- penguins[penguins$species == "Gentoo", ]
print(gentoo)
```
Tidyverse
``` r
gentoo  <- filter(penguins, species == "Gentoo")
print(gentoo)
# with pipe 
gentoo <- penguins %>%
    filter(species == "Gentoo")
print(gentoo)
```

Filtering and selecting in one with the pipe operator 
``` r
gentoo_bill_length <- penguins %>%
    filter(species == "Gentoo") %>%
    select(bill_length_mm)
print(gentoo_bill_length)
```
Now we understand a little bit of R we can move on to our first actual analysis. 

# Intro-to-RNAseq-analysis

---
title: "Intro-to-RNAseq-analysis"
author: "Adam Sorbie"
date: "2023-06-15"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Gut_brain_axis_course - Intro to RNAseq analysis

Today we will perform a fairly simple analysis of RNAseq data.This dataset comes
from another lab in our institute, AG Liesz, and was published a few years back: https://www.jneurosci.org/content/40/5/1162.long. Here the authors treated mice with short-chain fatty acids (produced by the gut microbiota) after stroke and found that it improved recovery.


Load the libraries: "DESeq2", "tidyverse", "EnhancedVolcano" and "clusterProfiler". 
Use the source function to source the "RNAseq_functions.R" file in the current working directory. 
```{r}
library(EnhancedVolcano)
library(DESeq2)
library(tidyverse)
library(clusterProfiler)
source("RNAseq_functions.R")
```
### Loading data 

Firstly we will read in the counts file and sample metadata using the ```read_tsv``` function from "readr", part of the "tidyverse". 
```{r}
counts <- read_tsv("GSE131788_counts_upload.SCFA_treatment.refseq_mm10.txt") %>% 
  # since we require our data to be numerical we need to move the first column to the rownames
  column_to_rownames("RefSeq") %>% 
  # deseq requires our data to be in matrix form 
  as.matrix()
metadata <- read_tsv("metadata_SCFA.txt")
```
Next we create a summarized experiment with the count data and sample metadata we just read in. A summarized experiment is an R object commonly used by bioconductor packages which simplifies and standardises storing genomic data in R.  
```{r}
se <- SummarizedExperiment(assays=list(counts=counts),colData = metadata)
```

In this analysis we will use the package DESeq2. DESeq2 is  one of the 
most common methods for detecting differentially expressed genes (DEGs) but also 
includes options for other types of analysis too. 
To load the data into DESeq2 to perform our analysis, we use the ```DESeqDataSet``` function, providing the summarized experiment and the design as a formula. In R we specify formulas using the ```~``` operator. Note that formulas do not require quotation. The variable on the left side is the dependent variable, while those on the right are the independent variables and are joined by ```+``` operators. 

e.g. Simple linear model - ```y ~ x```
     With two independent variables - ```y ~ x + b```
In our case the dependent variable or y, is the expression of a gene.     

Create a DESeq2 object called ```dds```using the DESeqDataSet function. This is
effectively 
```{r}
dds <- DESeqDataSet(se, design = ~ treatment)
```
Although not strictly necessary, we can also remove genes with low counts across samples to speed up computation and filter unimportant features. Since the assay is just a matrix we can extract it with the ```counts``` function and use ```rowSums``` to calculate which indices are above a given sum, in this case: 10. 
```{r}
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
```

### Principal component analysis

Often we want to see how similar or different our samples are. A simple way to
do this, is to use principal component analysis or PCA. PCA attempts to find a 
linear combination of variables which best explain the variation in the data.  
```{r}
vsd <- vst(dds, blind = T, fitType = "local")
plotPCA(vsd, intgroup=c("treatment")) + 
  theme_bw()
```

### Differential expression - which genes differ between conditions? 

Now we know our samples are different in terms of gene expression, we want
to find out which genes are altered between conditions.

To get our results we only need to run a single line of code. The ```DESeq```
function performs normalization, calculating 'size factors' to account for differences
in library depth and then estimate per-gene dispersion (variability). DESeq2 then 
uses a method called shrinkage to estimate more accurate estimates of variation which are
used in the final model. Finally, DESeq2 fits a negative binomial model to the data
and tests for differences in gene expression using a Wald test (default settings).

```{r}
dds <- DESeq(dds)
```
To access the DESeq results call the results function with the DESeq2 object as 
an argument.

```{r}
res <- results(dds)
```

Using the custom function ```annotate_degs``` let's add gene names and descriptions 
to our table. 

```{r}
res_annot <- annotate_degs(org.Mm.eg.db, res, keys = row.names(res), 
                           keytype = "REFSEQ", multivals = "first") %>% 
as.data.frame()
```

Finally, we can plot the results as a volcano plot. The "EnhancedVolcano" 
package makes this very easy but the defualt plots are not always the nicest/most 
interpretable

Plot volcano using default settings 
```{r}
EnhancedVolcano(res_annot,lab = res_annot$gene_symbol, x = "log2FoldChange", y = "padj")
```
Let's tidy it up a bit by adding a few extra arguments and generating key-value
pairs for the colour scheme so that SCFA enriched genes are coloured in red and 
those enriched in the control group are coloured in blue. 

Additionally, we can specify, foldchange and adjusted p-value thresholds, remove the
titles,captions and legend and then adjust the limits of the x-axis and now our
plot looks a lot better. 

```{r}
keyvals <- ifelse(
  res_annot$log2FoldChange < -1 & res_annot$padj < 0.1, 'dodgerblue',
    ifelse(res_annot$log2FoldChange > 1 & res_annot$padj < 0.1, 'firebrick1',
        'gray87'))
keyvals[is.na(keyvals)] <- 'gray87'
names(keyvals)[keyvals == 'firebrick1'] <- 'SCFA'
names(keyvals)[keyvals == 'dodgerblue'] <- 'Control'


EnhancedVolcano(res_annot, lab = res_annot$gene_symbol, x = "log2FoldChange", y = "padj", pCutoff = 0.1,
                FCcutoff = 1, xlim = c(-10, 10), title = NULL, subtitle = NULL, caption = NULL, 
                 colCustom = keyvals,colAlpha = 1, legendPosition = "None")
```
### Gene Ontology 

Although here we only have a handful of significantly regulated genes in many
cases we can have hundreds or even thousands. Even with properly annotated gene 
names, this can make functional interpretation difficult, thus we often assign 
functional categories to our gene lists to help make sense of them. For the last
section of todays course, we will perform overrepresentation analysis, assigning
biological profiles to our DEGs and using a hypergeometric test to determine
whether any of these categories are overrepresented compared to background. 

Firstly, we remove any NA values from our DESeq2 results and create a background
of all the genes included in our differential expression test. 
```{r}
res_annot_noNA <-  res_annot %>% 
  drop_na(padj)
## Create background dataset for hypergeometric testing using all tested genes for significance in the results                 
allOE_genes <- as.character(res_annot_noNA$gene_symbol)
```
Secondly, we extract the genes which were significantly different between
the SCFA and Control group. Since we have relatively few genes we are using a more
relaxed adjusted p-value filter of 0.1 and not filtering based on fold change. 
```{r}
## Extract significant results
sigOE_genes <- res_annot_noNA %>% 
  filter(padj < 0.1) %>% 
  pull(gene_symbol)
```
Next we will perform the enrichment analysis using the ```enrichGO``` function
from the "clusterProfiler" package. 
```{r}
## Run GO enrichment analysis 
ego <- enrichGO(gene = sigOE_genes, 
                universe = allOE_genes,
                keyType = "SYMBOL",
                OrgDb = org.Mm.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.1, 
                readable = TRUE)
```

Let's plot the results using the ```dotplot``` function, by providing our enrichment
results from above as an argument. 
```{r}
dotplot(ego)
```
