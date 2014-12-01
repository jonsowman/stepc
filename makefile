# A makefile for STEPC
#
# A fast dynamic system and controller simulator
# using Boost and uBLAS.
#
# Jon Sowman 2014 <jon@jonsowman.com>
#
# All Rights Reserved
#

GIT_VERSION := $(shell git describe --always --dirty)

CXX = g++
CXXFLAGS = -g -Wall -pedantic -DGIT_VERSION=\"$(GIT_VERSION)\"

TARGET = stepc

SOURCES  = $(wildcard *.cpp)
OBJS = $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

%.o: %.c
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY clean:
	rm -f $(TARGET) $(OBJS)
