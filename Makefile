SHELL = /bin/sh

INCLUDE_DIR := include
LIB_DIR := lib
SOURCE_DIR := src
BUILD_DIR := build

SOURCES	= $(wildcard $(SOURCE_DIR)/*.c)

CC	 = gcc
CXX = g++

CFLAGS	 = -I$(INCLUDE_DIR) -Wall -O3 -std=gnu99
CXXFLAGS = -fPIC -shared

LFLAGS	 = -lm 

ZLFLAGS = -L$(LIB_DIR) -linterface -lzcall -lstdc++ -lzkcm -lmpfr -lgmp -lgmpxx -fopenmp

all: libzcall.so libinterface.so backus_lsq backus_spr 

lib%.so: $(SOURCE_DIR)/%.cpp
	$(CXX) $< -o lib$* -I$(INCLUDE_DIR) $(ZLFLAGS) $(CXXFLAGS)

backus_lsq: bgv6_leastsq.c
	$(CC) -o $@ $^ $(CFLAGS) $(LFLAGS) $(ZLFLAGS)

backus_spr: bgv6_spread.c
	$(CC) -o $@ $^ $(CFLAGS) $(LFLAGS) $(ZLFLAGS)

