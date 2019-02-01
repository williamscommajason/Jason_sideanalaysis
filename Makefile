CC = g++
std=c++11

all: flacencode

flacencode: flacencode.o 
	$(CC) -std=$(std) flacencode.o -o flacencode -lFLAC

flacencode.o: flacencode.cpp
	$(CC) -std=$(std) -c flacencode.cpp


