#!/bin/bash

set -o errexit

CDIR="../../sode/cfiles"
OUTPUT="../cexamples.exe"

gcc -O3 -iquote "${CDIR}" -o "${OUTPUT}" \
           main.c \
           examples.c \
           "${CDIR}/randnorm.c" \
           "${CDIR}/solvers.c" 

echo "Compiled ${OUTPUT}"
