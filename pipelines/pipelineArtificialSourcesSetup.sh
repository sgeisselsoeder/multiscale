#!/bin/bash

# Make sure the binaries are up-to-date
make

### Prepare the events
inputFile=./input/ANTARES_2007-2012_random_generated_dummy/randomNeutrinoCandidatesInEquatorialCoordinates.txt
cp ${inputFile} ./input/events.txt

# Randomize all event positions
./bin/randomizeEvents --input ./input/events.txt
mv ./input/events.txtRandomized.txt ./input/events.txt
# Add artificial source events
./bin/addSourceEvents --input ./input/events.txt --dec 25.0 --ra 260.0 --radius 1.0 --number 10
mv ./input/events.txtWithSourceEvents.txt ./input/events.txt
./bin/addSourceEvents --input ./input/events.txt --dec 22.5 --ra 260.0 --radius 3.0 --number 30
mv ./input/events.txtWithSourceEvents.txt ./input/events.txt
./bin/addSourceEvents --input ./input/events.txt --dec 22.5 --ra 255.0 --radius 5.0 --number 50
mv ./input/events.txtWithSourceEvents.txt ./input/events.txt

