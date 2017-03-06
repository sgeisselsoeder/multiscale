#!/bin/bash

# Compute the expectation
time ./bin/step1 --input ./input/events.txt

# Evaluate the input distribution against the expectation
time ./bin/evaluateToSkymap --input ./input/events.txt
