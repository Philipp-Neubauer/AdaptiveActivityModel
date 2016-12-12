#!/bin/bash

set -ex

docker build -t docker.dragonfly.co.nz/squid-butterfish .

docker run --rm  -v $PWD:/work -w /work \
  docker.dragonfly.co.nz/squid-butterfish ./build.sh
