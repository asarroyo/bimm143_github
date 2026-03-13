# Class12
Areidy Arroyo (A17412951)

- [Section 1. Proportion of G/G in a
  population](#section-1-proportion-of-gg-in-a-population)
- [Section 4: Population Scale
  Analysis](#section-4-population-scale-analysis)

## Section 1. Proportion of G/G in a population

Download CSV file from Ensemble containing MXL from sample rs8067378

``` r
mxl <- read.csv("373531-SampleGenotypes-Homo_sapiens_Variation_Sample_rs8067378.csv")
head(mxl)
```

      Sample..Male.Female.Unknown. Genotype..forward.strand. Population.s. Father
    1                  NA19648 (F)                       A|A ALL, AMR, MXL      -
    2                  NA19649 (M)                       G|G ALL, AMR, MXL      -
    3                  NA19651 (F)                       A|A ALL, AMR, MXL      -
    4                  NA19652 (M)                       G|G ALL, AMR, MXL      -
    5                  NA19654 (F)                       G|G ALL, AMR, MXL      -
    6                  NA19655 (M)                       A|G ALL, AMR, MXL      -
      Mother
    1      -
    2      -
    3      -
    4      -
    5      -
    6      -

``` r
table(mxl$Genotype..forward.strand.)
```


    A|A A|G G|A G|G 
     22  21  12   9 

``` r
table(mxl$Genotype..forward.strand.) / nrow(mxl) * 100
```


        A|A     A|G     G|A     G|G 
    34.3750 32.8125 18.7500 14.0625 

To look for Great Britain we will look at the GBR population,still using
the CSV file we provided from earlier

``` r
gbr <- read.csv("373522-SampleGenotypes-Homo_sapiens_Variation_Sample_rs8067378.csv")
```

We will now find the proportion of G \| G

``` r
round(table(gbr$Genotype..forward.strand.) / nrow(gbr) * 100, 2)
```


      A|A   A|G   G|A   G|G 
    25.27 18.68 26.37 29.67 

## Section 4: Population Scale Analysis

We will be using multiple samples in order to be able to observe the
popultaion and its genetic differences based on population scale.

We have processed about 230 saomes and wille examine this to find its
relationship between rs8067378 and/or ORMDL3 expression

> Q13. Read this file into R and determine the sample size for each
> genotype and their corresponding median expression levels for each of
> these genotypes

Total samples

``` r
expr <- read.table("rs8067378_ENSG00000172057.6.txt")
head(expr)
```

       sample geno      exp
    1 HG00367  A/G 28.96038
    2 NA20768  A/G 20.24449
    3 HG00361  A/A 31.32628
    4 HG00135  A/A 34.11169
    5 NA18870  G/G 18.25141
    6 NA11993  A/A 32.89721

``` r
nrow(expr)
```

    [1] 462

``` r
table(expr$geno)
```


    A/A A/G G/G 
    108 233 121 

``` r
tapply(expr$exp, expr$geno, median)
```

         A/A      A/G      G/G 
    31.24847 25.06486 20.07363 

> Q14: Generate a boxplot with a box per genotype, what could you infer
> from the relative expression value between A/A and G/G displayed in
> this plot? Does the SNP effect the expression of ORMDL3?

Generate a box plot using `ggplot2()`

``` r
library(ggplot2)

ggplot(expr) + aes(x=geno, exp, fill=geno) +
geom_boxplot(notch=TRUE)
```

![](Class12_files/figure-commonmark/unnamed-chunk-10-1.png)

We can see that the SNP does affect ORMDL3 and its expression. We can
see that the genotypes A/A A/G and G/G all appear at different
expressions, some higher than others. Based on the boxplot we have
created the A/A genotype is much highly expressed than the others. We
can notice that G tends to decrease expression meanwhile A will increase
expression.
