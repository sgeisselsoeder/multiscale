#!/bin/bash

# Generate a random setup with artificial sources
./pipelineArtificialSourcesSetup.sh

# Compute the expectation
time ./bin/evaluateToSkymap --input ./input/events.txt

