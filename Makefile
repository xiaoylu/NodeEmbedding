CXX=g++
CXXFLAGS=-g -fopenmp -std=c++11 -Wall -pedantic 
BIN=sample

SRC=$(wildcard *.cpp)
OBJ=$(SRC:%.cpp=%.o)

%.o: %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

all: $(OBJ)
	$(CXX) -o $(BIN) $^ $(CXXFLAGS)

.PHONY: clean
clean:
	rm ./*.o

