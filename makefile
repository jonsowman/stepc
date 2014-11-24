# A makefile for STEPC
#
# A fast dynamic system and controller simulator
# using Boost and uBLAS.
#
# Jon Sowman 2014 <jon@jonsowman.com>
#
# All Rights Reserved
#
CXX = g++ -Wall -pedantic

TARGET = stepc

SOURCES  = $(wildcard *.cpp)
OBJS = $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS)

%.o: %.c
	$(CXX) -c $< -o $@

.PHONY clean:
	rm -f $(TARGET) $(OBJS)
