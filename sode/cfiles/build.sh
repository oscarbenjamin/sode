#!/bin/bash

set -o errexit

OUTPUT="../examples/cexamples.exe"

gcc -O3 -o "${OUTPUT}" \
           main.c \
           examples.c \
           "randnorm.c" \
           "solvers.c" \
           $@

echo "Compiled ${OUTPUT}"
