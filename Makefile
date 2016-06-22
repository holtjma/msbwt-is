CC=g++
CFLAGS=-O3 -std=c++11 -Wall
LDFLAGS=-O3 -std=c++11 -Wall
SOURCES=msbwtis.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=msbwtis

all: $(SOURCES) $(EXECUTABLE)
    
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) -c $(CFLAGS) $< -o $@
