#!/bin/bash

# Make sure the binaries are up-to-date
make

### Prepare the events
inputFile=./input/IceCube_2011-2014_random_generated_dummy/randomNeutrinoCandidatesInEquatorialCoordinates.txt
cp ${inputFile} ./input/events.txt

# Randomize all event positions
./bin/randomizeEvents --input ./input/events.txt
mv ./input/events.txtRandomized.txt ./input/events.txt
# Add artificial source events
./bin/addSourceEvents --input ./input/events.txt --dec 25.0 --ra 260.0 --radius 1.0 --number 30
cp ./input/events.txtWithSourceEvents.txt ./input/events.txt
./bin/addSourceEvents --input ./input/events.txt --dec 22.5 --ra 260.0 --radius 1.0 --number 30
cp ./input/events.txtWithSourceEvents.txt ./input/events.txt
./bin/addSourceEvents --input ./input/events.txt --dec 22.5 --ra 255.0 --radius 1.0 --number 30
# ./bin/addSourceEvents --input ./input/events.txt --dec 40 --ra 300 --decsize 4 --radius 8 --number 300
cp ./input/events.txtWithSourceEvents.txt ./input/events.txt

# Compute the expectation
date && qsub submitStep1
# wait until all previous computations have been finished
sleep 10 && while [ `qstat | wc -l` != 0 ] ; do qstat && sleep 900; done

# Produce pseudo experiments to assess relevances and significances later
date && for i in {1..100}; do qsub submitStep2 ; done
# wait until all previous computations have been finished
sleep 10 && while [ `qstat | wc -l` != 0 ] ; do qstat && sleep 900; done

date && time ./bin/step4
# date && qsub submitStep4
sleep 10 && while [ `qstat | wc -l` != 0 ] ; do qstat && sleep 900; done

date && qsub submitStep5

