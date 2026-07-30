CXX = g++
CXXFLAGS += -std=c++11 -Wall -pedantic
LDLIBS += -I$$mkEigenInc -fopenmp

EXE = main

all: build

.PHONY: all build clean distclean

build: $(EXE)

$(EXE): main.cpp
	$(CXX) $(CXXFLAGS) $(LDLIBS) $^ $(OUTPUT_OPTION)

clean:
	$(RM) *.o

distclean: clean
	$(RM) $(EXE)
	$(RM) *~