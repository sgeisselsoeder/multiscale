#!/bin/bash

# Make sure the binaries are up-to-date
make -C ..

### Prepare the events
inputFile=./input/IceCube_2011-2014_random_generated_dummy/randomNeutrinoCandidatesInEquatorialCoordinates.txt

echo "Evaluating data from $inputfile"
cp ../${inputFile} ../input/events.txt

echo "Compute the baseline expectation and wait until it is finished"
qsub submitStep1 
sleep 10 && while [ `qstat | wc -l` != 0 ] ; do qstat | wc -l && sleep 900; done

echo "Compute pseudo experiments to (in the end) deduce the significance"
for i in {1..100}; do qsub submitStep2 ; done
sleep 10 && while [ `qstat | wc -l` != 0 ] ; do qstat | wc -l && sleep 900; done

echo "Extract the distributions from the pseudo experiments"
qsub submitStep4
sleep 10 && while [ `qstat | wc -l` != 0 ] ; do qstat | wc -l && sleep 900; done

echo "Process the actual data, find clusters and deduce their significance compared to pseudo experiments"
qsub submitStep5
sleep 10 && while [ `qstat | wc -l` != 0 ] ; do qstat | wc -l && sleep 900; done

echo "Multiscale source search has finished. Thank you for your patience. You can find the results in the logfiles and in the subfolder intermediate/step5/."
