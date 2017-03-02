#!/bin/bash

here=`pwd`
# echo ${here}
for folder in `ls | grep scenario `; do echo ${folder} && cp gridLines.grid ${folder} && cp plotscript ${folder} && cd ${folder} && gnuplot plotscript ; cd ${here} ; done
