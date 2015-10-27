# This is configured for OS X to start. For Linux, HALIDE_CFLAGS probably needs
# to be extended with -lpthread -ldl
HALIDE_LIB ?= /usr/local/lib
HALIDE_INCLUDE ?= ${HALIDE_LIB}/../include
CXX ?= g++
HALIDE_CFLAGS ?= -std=c++11 -I ${HALIDE_INCLUDE} -L ${HALIDE_LIB} -lHalide

#need to fix this, not working
#-g is for debugging, check whether it should be in both statements
#%.standalone1: %.o
#	$(CXX) $(CFLAGS) -DSTANDALONE $< -g -o $@
#
#%.o: %.cpp %.h
#	$(CXX) $(CFLAGS) -DSTANDALONE $< -g -c $@

linearCombination_aot_compile: linearCombination_aot_compile.cpp
	$(CXX) $(HALIDE_CFLAGS) $< -g -o $@

lincombo_aot.o: linearCombination_aot_compile
	DYLD_LIBRARY_PATH=${HALIDE_LIB} ./linearCombination_aot_compile

linearCombination_aot_run: linearCombination_aot_run.cpp lincombo_aot.o
	$(CXX) $^ -g -o $@

run: linearCombination_aot_run
	./linearCombination_aot_run
