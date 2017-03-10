#!/bin/bash

set -ex

#Rscript -e 'install.packages("purrr",repo="http://cloud.r-project.org/")'
make 

cp draft.pdf /output/
