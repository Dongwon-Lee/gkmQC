CXX ?= g++
CXXFLAGS = -Wall -Wconversion -O3 -fPIC
LDFLAGS = -shared
SHVER = 2

SOURCES := $(shell find . -name '*.cpp')
HEADERS := $(shell find . -name '*.h')
OBJECTS = $(patsubst $.cpp,obj/%.o$(SOURCES))
OUTPUT  := GkmKernel.so

OBJ_GKMLIB = GkmKernel.o

all: $(OUTPUT)
$(OUTPUT): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $(SOURCES) -lm -lpthread $(LDFLAGS) -o $(OUTPUT) 

clean:
	rm $(OUTPUT)