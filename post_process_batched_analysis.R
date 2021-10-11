#!/usr/bin/Rscript --vanilla

# post-process batched analysis

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

filestoload = args[[1]]

data = list()
for(f in filestoload){
	data = readRDS()
}