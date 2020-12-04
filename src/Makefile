CXX ?= g++
CFLAGS = -Wall -Wconversion -O3 -fPIC
LDFLAGS = -shared
SHVER = 2

all: gkmkern 

gkmkern: gkmkern_main.c libgkm.o
	$(CXX) $(CFLAGS) gkmkern_main.c libgkm.o -o gkmkern -lm -lpthread

libgkm.so: libgkm.o
	$(CXX) $(CFLAGS) libgkm.o -lm -lpthread -o $@ $(LDFLAGS)

libgkm.o: libgkm.c libgkm.h
	$(CXX) $(CFLAGS) -c libgkm.c

clean:
	rm -f libgkm.o gkmkern

install:
	cp libgkm.so gkmkern ../bin
