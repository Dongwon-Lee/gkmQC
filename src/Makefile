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

# install:
#   - ``gkmkern_pylib.so`` goes next to the Python package so
#     ``gkmsvm.py`` can ``ctypeslib.load_library`` it via
#     ``os.path.dirname(__file__)`` regardless of how the package was
#     later installed (editable / wheel / source tree).
#   - ``gkmkern`` (the standalone CLI) goes on PATH: ``$CONDA_PREFIX/bin``
#     when running under an active conda env, else ``../bin`` — the
#     repo-local fallback that kept the pre-PR README flow working.
install:
	cp gkmkern_pylib.so ../gkmqc/
	@if [ -n "$$CONDA_PREFIX" ] && [ -d "$$CONDA_PREFIX/bin" ]; then \
	    echo "installing gkmkern into $$CONDA_PREFIX/bin"; \
	    cp gkmkern $$CONDA_PREFIX/bin/; \
	else \
	    echo "installing gkmkern into ../bin (no CONDA_PREFIX detected)"; \
	    cp gkmkern ../bin/; \
	fi
