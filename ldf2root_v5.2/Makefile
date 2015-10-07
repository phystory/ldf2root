
all : ldf2root

LIBS:=$(shell root-config --libs)
INCDIR:=$(shell root-config --incdir)
GCC=g++

ldf2root :  ldf2root.o 
	$(GCC) -O ldf2root.o $(LIBS)  -o ldf2root

ldf2root.o :  ldf2root.cpp 
	$(GCC) -g -O -Wall -fPIC -pthread -std=c++11 -I$(INCDIR) -c ldf2root.cpp

clean:
	rm -f *o
	rm -f ldf2root
	rm -f readroot.C
