TARGET = ising

HDF5_ROOT=/usr/local
HDF5_INC=$(HDF5_ROOT)/include
HDF5_LIB=$(HDF5_ROOT)/lib


CC = gcc-9
CPPFLAGS += -I./include -I$(HDF5_INC)
CFLAGS= -fopenmp
LIBS = -lm $(HDF5_LIB) -lhdf5

.PHONY: clean all default

default: $(TARGET)
all: default

OBJECTS = src/ising_functions.o src/ising.o
HEADERS = $(wildcard *.h)

ising_functions.o: src/ising_functions.c
	$(CC) $(CPPFLAGS)  $(CFLAGS) -c $< -o $@

ising.o: src/ising.c
	$(CC) $(CPPFLAGS) $(CFLAGS)  -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OBJECTS) -Wall $(LIBS) -o $@

clean:
	-rm -f src/*.o
	-rm -f $(TARGET)
