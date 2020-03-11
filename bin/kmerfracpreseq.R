#!/usr/bin/env Rscript --vanilla
# A script to estimate genome coverage from a tabular, two column kmer histogram file  
# Adam Rivers, USDA-ARS-GBRU,  July 13, 2017

## Load Libraries ##
# Modified version of https://github.com/smithlabcode/preseqR ver. 3.1.1 with kmer.frac.bootstrap, kmer.frac.curve exported to Namespase 
if(!require("preseqR")){
  install.packages("../assets/preseqR.tar.gz",  repo=NULL, type="source")
  library("preseqR")
}
library(preseqR)
library(argparser, quietly=TRUE)

# Create a parser
p <- arg_parser("A script to estimate genome coverage from a tabular, two column kmer histogram file, creating a coverage table")
# Add command line arguments)
p <- add_argument(p, "--input", help="Two column histogram input file", type="character")
p <- add_argument(p, "--output", help="Output file path", type="character")
p <- add_argument(p, "--points", help="The number of points to evaluate coverage over", default=100)
p <- add_argument(p, "--bases", help="The number of bases in the sample, if omitted, calculated from kmer data",type="numeric", default=NA)
p <- add_argument(p, "--coverage", help="coverage depths to evaluate, a vector of positive integers e.g. c(2, 5, 10)", type="numeric", default=c(2,5,10))


# Parse the command line arguments 
args <- parse_args(p)


## Functions ##

# create log spaced points
logspace <- function( d1, d2, n){
  exp(log(10)*seq(log10(d1), log10(d2), length.out=n)) 
}

# create a list of log spaced points around the sequencing depth of the sample (in bases)
create.size<- function(bases,points){
  min=bases/10000
  points=logspace(min,bases*100,points)
  return(points)
}

## Main ##
# load input kmer histogram data
n <- read.table(args$input, sep="\t")

# if no --bases value was supplied calulate from the data
if (is.na(args$bases)){
  args$bases <- sum(n[,2])
  warning("Value not supplied for --bases, using the sum of the Kmers from input file")
}

#  create a vector of coverage points to evaluate over in GB
size <- create.size(bases=args$bases,points=args$points)/1e9

# Run the rational function approxamation 
curves <- kmer.frac.curve(n=n, N=args$bases, seq.size.GB=size, r=args$coverage, mt=100)

#write the results
write.table(x=curves, file=args$output, row.names=FALSE)