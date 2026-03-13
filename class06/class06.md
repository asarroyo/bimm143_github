# Class06: R Functions
Areidy (PID: A17412951)

- [Background](#background)
- [A first function](#a-first-function)
- [A second function](#a-second-function)

## Background

Functions are at the heart of using R. Everything we do involves calling
and using function (from data input, analysis to results output)

All functions in R have at least 3 things:

1.  A **name** the thing we use to call the function.
2.  One or more input **arguments** that a are comma separated
3.  The **body**, lines of code between curly brackets { } that does the
    work of the function.

## A first function

Let’s write a silly wee function to add some number:

``` r
add <- function(x) {
  x + 1
}
```

Let’s try it out

``` r
add(100)
```

    [1] 101

Will this work

``` r
add(c(100, 200, 300))
```

    [1] 101 201 301

Modify to be more useful and add more than just 1

``` r
add <- function(x, y=1) {
  x + y
}
```

``` r
add(100, 10)
```

    [1] 110

Will this still work

``` r
add(100)
```

    [1] 101

``` r
plot(1:10, col="blue", typ="b")
```

![](class06_files/figure-commonmark/unnamed-chunk-7-1.png)

``` r
log(10, base=10)
```

    [1] 1

> **N.B.** Input arguments can be either **required** or **optional**.
> The later have a fall-back default that is specified in the function
> code with an equals sign.

``` r
#add(x=100, y=200, z=300)
```

## A second function

All functions in R look like this

    name <- function(arg) { 
      body 
    }

The `sample()` function in R …

``` r
sample(1:10, size = 4)
```

    [1] 6 8 2 3

> Q. Return 12 numbers picked randomly from the input 1:10

``` r
sample(1:10, size = 12, replace = TRUE)
```

     [1] 9 3 8 9 5 6 3 7 3 1 8 7

> Q. Write the code to generate a random 12 nucelotide long DNA
> sequence?

``` r
bases <- c("A", "C", "G", "T")
sample(bases, size=12, replace=TRUE)
```

     [1] "T" "C" "T" "G" "C" "T" "G" "T" "G" "G" "C" "C"

> Q. Write a first version function called `generate_dna()` that
> generates a user specified length `n` random DNA sequence?

``` r
generate_dna <- function(n=6) { 
  bases <- c("A", "C", "G", "T")
  sample(bases, size = n, replace = TRUE)
}
```

``` r
generate_dna(100)
```

      [1] "A" "T" "A" "A" "G" "G" "A" "C" "G" "A" "A" "C" "A" "T" "A" "G" "T" "T"
     [19] "A" "T" "G" "G" "T" "C" "C" "G" "A" "A" "C" "T" "C" "A" "A" "C" "G" "T"
     [37] "G" "C" "G" "A" "C" "T" "C" "T" "A" "A" "G" "G" "G" "T" "C" "T" "G" "T"
     [55] "A" "A" "C" "G" "G" "G" "T" "T" "C" "T" "A" "T" "C" "A" "G" "A" "T" "T"
     [73] "T" "T" "G" "C" "A" "A" "T" "C" "G" "A" "G" "G" "T" "G" "G" "G" "T" "A"
     [91] "C" "A" "C" "T" "A" "G" "T" "A" "C" "G"

> Q. Modify your funtion to return a FASTA like sequence so rather than
> \[1\] “T” “C” “T” “G” “C” we want “TCTGC”

``` r
generate_dna <- function(n=6) { 
  bases <- c("A", "C", "G", "T")
  ans <- sample(bases, size = n, replace = TRUE)
  ans <- paste(ans, collapse = "")
  return(ans)
  x <- "pooopoppants"
}
```

``` r
generate_dna(10)
```

    [1] "CGGCGGGGTT"

An example

``` r
# Example pattern (not using your bases)
x <- c("T", "C", "T", "G", "C")

paste(x, collapse = "TCTGC")
```

    [1] "TTCTGCCTCTGCTTCTGCGTCTGCC"

``` r
# returns "HELLO"
```

> . Q. Give the user an option to return to FASTA format output sequence
> or standard multi-elemnt vector format?

``` r
generate_dna <- function(n=6, fasta= TRUE) { 
  bases <- c("A", "C", "G", "T")
  ans <- sample(bases, size = n, replace = TRUE)
  
  if (fasta) { 
    ans <- paste(ans, collapse = "") 
    cat("Hello...")
  } else {
      cat("...is it me you are looking for...")
  }

  return(ans)
}
```

``` r
generate_dna(10)
```

    Hello...

    [1] "GCGCATGTAA"

``` r
generate_dna(10, fasta = F)
```

    ...is it me you are looking for...

     [1] "A" "C" "G" "C" "A" "A" "C" "A" "C" "C"

> Q. Write a function called `generate_protein()` that generates a user
> specified length protein sequence in FASTA like format?

``` r
generate_protein <- function(n){
  
aa <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  
  ans <- sample(aa, size = n, replace = T)
  ans <- paste (ans, collapse= "")
  return(ans)
}
```

``` r
generate_protein(10)
```

    [1] "YFGGAFWFLR"

> Q. Use your new `generate_protein()` function to generate all
> sequences length 6 and 12 amin-acids in length and check of any of
> these are unique in nature(i.e. found in the NR database at NCBI)?

``` r
generate_protein(6)
```

    [1] "WRMLGN"

``` r
generate_protein(7)
```

    [1] "SCGWAHF"

``` r
generate_protein(8)
```

    [1] "HERPHLGM"

``` r
generate_protein(9)
```

    [1] "ADKTIEYIY"

``` r
generate_protein(10)
```

    [1] "RMWQSTLIPP"

``` r
generate_protein(11)
```

    [1] "KWYMDYPPPWL"

``` r
generate_protein(12)
```

    [1] "PKWQKHWWAWIH"

or we could do a `for()` loop:

``` r
for(i in 6:12) {
  cat(">", i, sep= "", "\n")
  cat( generate_protein(i), "\n")
}
```

    >6
    WYAQDE 
    >7
    YRPLWFS 
    >8
    TENREFSF 
    >9
    QCQQMSIFG 
    >10
    RYGGSMSFFR 
    >11
    HPFTKCSVANK 
    >12
    WVTWTCQLAMNI 

> id AGKRTST next G

N.B\> 6-8 sequences were found on NCBI
