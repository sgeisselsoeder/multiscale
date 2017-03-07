#!/bin/bash

# create/update the binaries
make

# Generate a random setup with artificial sources
./pipelines/setupArtificialSources.sh

# Compute the expectation
time ./bin/step1 --input ./input/events.txt

# Evaluate the input distribution against the expectation
time ./bin/evaluateToSkymap --input ./input/events.txt

cd ./output/
cp plotscript scenarioEvaluation/
cp gridLines.grid scenarioEvaluation/
cd scenarioEvaluation/
gnuplot plotscript

