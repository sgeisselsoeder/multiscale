#!/bin/bash

echo "#### Computing the multiscale search on this machine locally. This is going to take some time ..."
echo "#### Please make sure that you execute this pipeline scripts from within the pipeline folder, \
	namely via ./pipelineMultiscaleLocal.sh, as the relative paths are set accordingly."
cd ..

echo "Make sure the binaries are up-to-date"
make

inputFile=./input/IceCube_2011-2014_random_generated_dummy/randomNeutrinoCandidatesInEquatorialCoordinates.txt
echo "#### Evaluating data from" ${inputFile}
cp ./${inputFile} ./input/events.txt

echo "#### Compute the baseline expectation"
time ./bin/step1 --input "./input/events.txt"

echo "#### Compute pseudo experiments to (in the end) deduce the significance"
for i in {1..10}; do echo "#### Iteration " $i " of 10" && ./bin/step2 --input "./input/events.txt" ; done

echo "#### Extract the distributions from the pseudo experiments"
time ./bin/step4

echo "#### Process the actual data, find clusters and deduce their significance compared to pseudo experiments"
time ./bin/step5 --input "./input/events.txt"

echo "#### Multiscale source search has finished. Thank you for your patience. You can find the results in the logfiles and in the subfolder intermediate/step5/."
