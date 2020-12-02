CXX ?= g++
CFLAGS = -Wall -Wconversion -O3 -fPIC
#CFLAGS = -Wall -Wconversion -g -fPIC
LDFLAGS = -shared
SHVER = 2

all: GkmKernel.so

GkmKernel.so: libgkm.o GkmKernel.o
	$(CXX) $(CFLAGS) libgkm.o GkmKernel.o -lm -lpthread -o $@ $(LDFLAGS)

#libsvm.o: libsvm.cpp libsvm.h
#	$(CXX) $(CFLAGS) -c libsvm.cpp
GkmKernel.o: GkmKernel.c GkmKernel.h
	$(CXX) $(CFLAGS) -c GKMKernel.c

clean:
	rm -f GkmKernel.o GkmKernel.so 

install:
	cp GkmKernel.so ../bin