#!/bin/bash

# this pipeline tests splitting available data to derive a hypothesis on one part of the data and trying to confirm it in the second part

# Make sure the binaries are up-to-date
make

### Prepare the events
inputFile=./input/IceCube_2011-2014_random_generated_dummy/randomNeutrinoCandidatesInEquatorialCoordinates.txt
cp ${inputFile} ./input/events.txt
./bin/randomizeEvents --input ./input/events.txt
mv ./input/events.txtRandomized.txt ./input/events.txt
cp ./input/events.txt ./input/eventsBackup.txt

# Add artificial source events
./bin/addSourceEvents --input ./input/events.txt --dec 70.0 --ra 35.0 --radius 10.0 --decsize 10.0 --number 400
cp ./input/events.txtWithSourceEvents.txt ./input/events.txt

# Prepare two data subsets, one for the generation of hypotheses and one for the test
./bin/shuffleEvents --input ./input/events.txt
cp ./input/events.txtShuffled.txt ./input/events.txt
./bin/splitEvents --input ./input/events.txt --splitRatio 0.5

### Evaluate the first data subset to generate hypotheses
cp ./input/events.txtSplitPart1.txt ./input/events.txt
# Compute the expectation and wait until it is finished
qsub submitStep1 && sleep 10 && while [ `qstat | wc -l` != 0 ] ; do qstat && sleep 900; done
qsub submitEvaluation
# Wait for all jobs to finish
while [ `qstat | wc -l` != 0 ] ; do qstat && sleep 900; done
# Save the in and outputs for later analysis
cp -r output/scenarioEvaluation output/scenarioEvaluationHypothesisGeneration
cp -r intermediate/step1/normalizationDataStep1/ intermediate/step1/normalizationDataStep1HypothesisGeneration/

### Evaluate the second data subset to confirm the hypotheses
cp ./input/events.txtSplitPart2.txt ./input/events.txt
# Compute the expectation and wait until it is finished
qsub submitStep1 && sleep 10 && while [ `qstat | wc -l` != 0 ] ; do qstat && sleep 900; done 
qsub submitEvaluation
# Wait for all jobs to finish
while [ `qstat | wc -l` != 0 ] ; do qstat && sleep 900; done
# Save the in and outputs for later analysis
cp -r output/scenarioEvaluation output/scenarioEvaluationHypothesisTesting
cp -r intermediate/step1/normalizationDataStep1/ intermediate/step1/normalizationDataStep1HypothesisTesting

# Produce statistics from pseudo experiments to be able to assess the significance(s)
for i in {1..50}; do qsub submitEvaluationRandom ; done



