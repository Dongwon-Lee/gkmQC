CXX ?= g++
CFLAGS = -Wall -Wconversion -O3 -fPIC
LDFLAGS = -shared
SHVER = 2

all: gkmkern gkmkern_pylib.so

# making executive binary
gkmkern: gkmkern_main.c libgkm.o
	$(CXX) $(CFLAGS) gkmkern_main.c libgkm.o -o gkmkern -lm -lpthread

# making shared binary
gkmkern_pylib.so: gkmkern_pylib.o libgkm.o
	$(CXX) $(CFLAGS) gkmkern_pylib.o libgkm.o -lm -lpthread -o $@ $(LDFLAGS)

libgkm.o: libgkm.c libgkm.h
	$(CXX) $(CFLAGS) -c libgkm.c

gkmkern_pylib.o: gkmkern_pylib.c 
	$(CXX) $(CFLAGS) -c gkmkern_pylib.c

clean:
	rm -f libgkm.o gkmkern_pylib.o gkmkern_pylib.so gkmkern

install:
	cp gkmkern_pylib.so gkmkern ../bin
