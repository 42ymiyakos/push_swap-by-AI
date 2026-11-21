# Research build. `make` builds the solver (psig), the SA prototype (qsa) and the checker.
CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra
all: bin/psig bin/qsa tools/checker
bin/psig: src/psig.cpp src/*.hpp
	@mkdir -p bin; $(CXX) $(CXXFLAGS) -o $@ src/psig.cpp
bin/qsa: src/qsa.cpp src/*.hpp
	@mkdir -p bin; $(CXX) $(CXXFLAGS) -o $@ src/qsa.cpp
tools/checker: tools/checker.c
	$(CC) -O2 -Wall -Wextra -o $@ $<
asan: src/psig.cpp src/*.hpp
	@mkdir -p bin; $(CXX) -O1 -g -std=c++17 -fsanitize=address,undefined -o bin/psig_asan src/psig.cpp
sets:
	tools/make_sets.sh
clean:
	rm -f bin/psig bin/qsa bin/psig_asan tools/checker
.PHONY: all asan sets clean
