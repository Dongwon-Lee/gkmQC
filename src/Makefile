CXX ?= g++
CFLAGS = -Wall -Wconversion -O3 -fPIC
LDFLAGS = -shared
SHVER = 2

all: libgkm.so

libgkm.so: libgkm.o
	$(CXX) $(CFLAGS) libgkm.o -lm -lpthread -o $@ $(LDFLAGS)

libgkm.o: libgkm.c libgkm.h
	$(CXX) $(CFLAGS) -c libgkm.c

clean:
	rm -f libgkm.o libgkm.so 

install:
	cp libgkm.so ../bin