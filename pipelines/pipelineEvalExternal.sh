#!/bin/bash

# This pipeline makes no sense if you don't have a certain follow-up project ...

# Make sure the binaries are up-to-date
make -C ..

### Prepare the events
inputFile=./input/nonpublicData/ICFull/IceCubeDatasetLocalCoords.txt
cp ../${inputFile} ../input/events.txt

# Compute the expectation and wait until it is finished
qsub submitStep1 
sleep 10 && while [ `qstat | wc -l` != 0 ] ; do qstat && sleep 900; done

# Evaluate the actual data to skymap
qsub submitEvaluation

# Produce statistics from pseudo experiments to be able to assess the significance(s)
for i in {1..100}; do qsub submitEvaluationRandom ; done

