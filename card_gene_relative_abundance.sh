#!/bin/bash

# compute relative abundances of card genes

rgi bwt --include_wildcard -1 $2 -2 $3 -a bwa -n 10 -o "$1"_rgibwt/"$1"_rgibwt  
