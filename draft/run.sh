#!/bin/bash

set -ex

docker build -t docker.dragonfly.co.nz/auckland-bivales .

docker run --rm  -v $PWD:/work -w /work \
  docker.dragonfly.co.nz/auckland-bivales ./build.sh
