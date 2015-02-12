
all : ldf2root

LIBS=-L$(ROOTLIB) -lCore  -lCint -lHist -lGraf -lGraf3d -lGpad -lTree -lRIO -lRint -lPostscript -lMatrix -lPhysics -pthread -lm -ldl -rdynamic
GCC=g++

ldf2root :  ldf2root.o 
	$(GCC) -O ldf2root.o $(LIBS)  -o ldf2root

ldf2root.o :  ldf2root.cpp 
	$(GCC) -g -O -Wall -fPIC -pthread -I$(ROOTINCLUDE) -c ldf2root.cpp

clean:
	rm -f *o
	rm -f ldf2root
	rm -f readroot.C
