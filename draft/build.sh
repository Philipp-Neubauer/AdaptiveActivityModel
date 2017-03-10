#!/bin/bash

set -ex

#Rscript -e 'install.packages("purrr",repo="http://cloud.r-project.org/")'
make -B

cp draft.pdf /output/
