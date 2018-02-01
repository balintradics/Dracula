#!/bin/bash

g++ `gsl-config --cflags` -c Dracula.cxx
g++ `gsl-config --cflags` -c main.cxx
g++ -o main main.o Dracula.o `gsl-config --libs` -lsqlite3
#g++ `gsl-config --cflags` -c mainFlight.cxx
#g++ -o mainFlight mainFlight.o Dracula.o `gsl-config --libs` -lsqlite3
