#!/bin/bash

set -ex

docker build -t docker.dragonfly.co.nz/adaptive_activity .

docker run --rm  -v $PWD:/work -w /work \
  docker.dragonfly.co.nz/adaptive_activity ./build.sh
