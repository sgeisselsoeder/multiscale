CXX      = g++
CXXFLAGS  = -std=c++11 -Wall -Wextra -pedantic -fopenmp -Wunreachable-code # -Werror # -Weffc++
INCFLAGS = -Iinclude -I${BOOST_INCDIR}/
LDFLAGS = -L${BOOST_LIBDIR}/ -lboost_program_options

CXXFLAGS += ${INCFLAGS}

### Debug
# CXXFLAGS += -g -fprofile-arcs -ftest-coverage
# LDFLAGS += -fprofile-arcs -lgcov
### Release
CXXFLAGS += -O3 -DNDEBUG
 

SRCS   = $(wildcard src/*.cpp)
TARGETTESTS   = tests
TARGETSTEP1   = step1
TARGETSTEP2   = step2
TARGETSTEP4   = step4
TARGETSTEP5   = step5
TARGETSTEP5RANDOM   = evaluateRandom
TARGETSTEP5RANDOMWSOURCE   = evaluateRandomWithSource
TARGETRANDOMIZEEVENTS    = randomizeEvents
TARGETADDSOURCEEVENTS    = addSourceEvents
TARGETSPLITEVENTS    = splitEvents
TARGETSHUFFLEEVENTS    = shuffleEvents
TARGETTOGETHER   = together
TARGETEVAL   = evaluateToSkymap
OBJS   = $(SRCS:.cpp=.o)
DEPS   = $(SRCS:.cpp=.depends)

.PHONY: clean all doc purge

all: $(TARGETTESTS) $(TARGETSTEP1) $(TARGETSTEP2) $(TARGETSTEP4) $(TARGETSTEP5) $(TARGETSTEP5RANDOM) $(TARGETSTEP5RANDOMWSOURCE) $(TARGETTOGETHER) $(TARGETEVAL) $(TARGETADDSOURCEEVENTS) $(TARGETRANDOMIZEEVENTS) $(TARGETSPLITEVENTS) $(TARGETSHUFFLEEVENTS) 

clean:
	rm -f $(OBJS) $(DEPS) src/bins/*.o bin/*
	rm -rf doc/html doc/latex
	find . -iname "*.gc*" -exec rm -rf {} \;
	
purge: clean
	rm -f intermediate/step*/*.*
	rm -f intermediate/step*/*/*.*
	# rm -f output/*.*
	rm -f output/*/*
	
re: clean all
	echo "Refreshed"
	
$(TARGETSHUFFLEEVENTS): $(OBJS) src/bins/$(TARGETSHUFFLEEVENTS).o
	$(CXX) $(CXXFLAGS) $(OBJS) src/bins/$(TARGETSHUFFLEEVENTS).o -o bin/$(TARGETSHUFFLEEVENTS) $(LDFLAGS)

$(TARGETSHUFFLEEVENTS).o: src/bins/$(TARGETSHUFFLEEVENTS).cpp
	$(CXX) $(CXXFLAGS) -c src/bins/$(TARGETSHUFFLEEVENTS).cpp -o src/bins/$(TARGETSHUFFLEEVENTS).o

$(TARGETSPLITEVENTS): $(OBJS) src/bins/$(TARGETSPLITEVENTS).o
	$(CXX) $(CXXFLAGS) $(OBJS) src/bins/$(TARGETSPLITEVENTS).o -o bin/$(TARGETSPLITEVENTS) $(LDFLAGS)

$(TARGETSPLITEVENTS).o: src/bins/$(TARGETSPLITEVENTS).cpp
	$(CXX) $(CXXFLAGS) -c src/bins/$(TARGETSPLITEVENTS).cpp -o src/bins/$(TARGETSPLITEVENTS).o

$(TARGETRANDOMIZEEVENTS): $(OBJS) src/bins/$(TARGETRANDOMIZEEVENTS).o
	$(CXX) $(CXXFLAGS) $(OBJS) src/bins/$(TARGETRANDOMIZEEVENTS).o -o bin/$(TARGETRANDOMIZEEVENTS) $(LDFLAGS)

$(TARGETRANDOMIZEEVENTS).o: src/bins/$(TARGETRANDOMIZEEVENTS).cpp
	$(CXX) $(CXXFLAGS) -c src/bins/$(TARGETRANDOMIZEEVENTS).cpp -o src/bins/$(TARGETRANDOMIZEEVENTS).o
	
$(TARGETADDSOURCEEVENTS): $(OBJS) src/bins/$(TARGETADDSOURCEEVENTS).o
	$(CXX) $(CXXFLAGS) $(OBJS) src/bins/$(TARGETADDSOURCEEVENTS).o -o bin/$(TARGETADDSOURCEEVENTS) $(LDFLAGS)

$(TARGETADDSOURCEEVENTS).o: src/bins/$(TARGETADDSOURCEEVENTS).cpp
	$(CXX) $(CXXFLAGS) -c src/bins/$(TARGETADDSOURCEEVENTS).cpp -o src/bins/$(TARGETADDSOURCEEVENTS).o

$(TARGETEVAL): $(OBJS) src/bins/$(TARGETEVAL).o
	$(CXX) $(CXXFLAGS) $(OBJS) src/bins/$(TARGETEVAL).o -o bin/$(TARGETEVAL) $(LDFLAGS)

$(TARGETEVAL).o: src/bins/$(TARGETEVAL).cpp
	$(CXX) $(CXXFLAGS) -c src/bins/$(TARGETEVAL).cpp -o src/bins/$(TARGETEVAL).o

$(TARGETTESTS): $(OBJS) src/bins/$(TARGETTESTS).o
	$(CXX) $(CXXFLAGS) $(OBJS) src/bins/$(TARGETTESTS).o -o bin/$(TARGETTESTS) $(LDFLAGS)

$(TARGETTESTS).o: src/bins/$(TARGETTESTS).cpp
	$(CXX) $(CXXFLAGS) -c src/bins/$(TARGETTESTS).cpp -o src/bins/$(TARGETTESTS).o
	
$(TARGETTOGETHER): $(OBJS) src/bins/$(TARGETTOGETHER).o
	$(CXX) $(CXXFLAGS) $(OBJS) src/bins/$(TARGETTOGETHER).o -o bin/$(TARGETTOGETHER) $(LDFLAGS)

$(TARGETTOGETHER).o: src/bins/$(TARGETTOGETHER).cpp
	$(CXX) $(CXXFLAGS) -c src/bins/$(TARGETTOGETHER).cpp -o src/bins/$(TARGETTOGETHER).o

$(TARGETSTEP1): $(OBJS) src/bins/$(TARGETSTEP1).o
	$(CXX) $(CXXFLAGS) $(OBJS) src/bins/$(TARGETSTEP1).o -o bin/$(TARGETSTEP1) $(LDFLAGS)

$(TARGETSTEP1).o: src/bins/$(TARGETSTEP1).cpp
	$(CXX) $(CXXFLAGS) -c src/bins/$(TARGETSTEP1).cpp -o src/bins/$(TARGETSTEP1).o

$(TARGETSTEP2): $(OBJS) src/bins/$(TARGETSTEP2).o
	$(CXX) $(CXXFLAGS) $(OBJS) src/bins/$(TARGETSTEP2).o -o bin/$(TARGETSTEP2) $(LDFLAGS)

$(TARGETSTEP2).o: src/bins/$(TARGETSTEP2).cpp
	$(CXX) $(CXXFLAGS) -c src/bins/$(TARGETSTEP2).cpp -o src/bins/$(TARGETSTEP2).o
	
$(TARGETSTEP4): $(OBJS) src/bins/$(TARGETSTEP4).o
	$(CXX) $(CXXFLAGS) $(OBJS) src/bins/$(TARGETSTEP4).o -o bin/$(TARGETSTEP4) $(LDFLAGS)

$(TARGETSTEP4).o: src/bins/$(TARGETSTEP4).cpp
	$(CXX) $(CXXFLAGS) -c src/bins/$(TARGETSTEP4).cpp -o src/bins/$(TARGETSTEP4).o
	
$(TARGETSTEP5): $(OBJS) src/bins/$(TARGETSTEP5).o
	$(CXX) $(CXXFLAGS) $(OBJS) src/bins/$(TARGETSTEP5).o -o bin/$(TARGETSTEP5) $(LDFLAGS)

$(TARGETSTEP5).o: src/bins/$(TARGETSTEP5).cpp
	$(CXX) $(CXXFLAGS) -c src/bins/$(TARGETSTEP5).cpp -o src/bins/$(TARGETSTEP5).o

$(TARGETSTEP5RANDOM): $(OBJS) src/bins/$(TARGETSTEP5RANDOM).o
	$(CXX) $(CXXFLAGS) $(OBJS) src/bins/$(TARGETSTEP5RANDOM).o -o bin/$(TARGETSTEP5RANDOM) $(LDFLAGS)

$(TARGETSTEP5RANDOM).o: src/bins/$(TARGETSTEP5RANDOM).cpp
	$(CXX) $(CXXFLAGS) -c src/bins/$(TARGETSTEP5RANDOM).cpp -o src/bins/$(TARGETSTEP5RANDOM).o

$(TARGETSTEP5RANDOMWSOURCE): $(OBJS) src/bins/$(TARGETSTEP5RANDOMWSOURCE).o
	$(CXX) $(CXXFLAGS) $(OBJS) src/bins/$(TARGETSTEP5RANDOMWSOURCE).o -o bin/$(TARGETSTEP5RANDOMWSOURCE) $(LDFLAGS)

$(TARGETSTEP5RANDOMWSOURCE).o: src/bins/$(TARGETSTEP5RANDOMWSOURCE).cpp
	$(CXX) $(CXXFLAGS) -c src/bins/$(TARGETSTEP5RANDOMWSOURCE).cpp -o src/bins/$(TARGETSTEP5RANDOMWSOURCE).o

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.depends: %.cpp
	$(CXX) -M $(CXXFLAGS) $< > $@
	
doc:
	doxygen doxygen.config
	
-include $(DEPS)

