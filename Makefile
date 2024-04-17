SHELL := bash

INCLUDE_DIR := include
LIB_DIR := lib
SOURCE_DIR := src

SOURCES=$(wildcard $(SOURCE_DIR)/*.c)
UTILS=$(SOURCES:$(SOURCE_DIR)/%.c=%)

CC  = gcc
CXX = g++

CFLAGS	 = -I$(INCLUDE_DIR) -L$(LIB_DIR) -Wall  -std=gnu99
CXXFLAGS = -I$(INCLUDE_DIR) -L$(LIB_DIR) -fPIC -shared

LFLAGS	 = -lm 

ZLFLAGS = -lstdc++ -lzkcm -lmpfr -lgmp -lgmpxx -fopenmp

all: libzcall.so libinterface.so backus_hmr $(UTILS)

# Compile libraries

libzcall.so: $(SOURCE_DIR)/zcall.cpp
	mkdir -p $(LIB_DIR)
	$(CXX) $< $(ZLFLAGS) -o $(LIB_DIR)/$@ $(CXXFLAGS)

libinterface.so: $(SOURCE_DIR)/interface.cpp
	mkdir -p $(LIB_DIR)
	$(CXX) $< $(ZLFLAGS) -o $(LIB_DIR)/$@ $(CXXFLAGS) -lzcall

# Compile utility scripts

%: $(SOURCE_DIR)/%.c
	$(CC) $^ -o $@ $(LFLAGS)

# Compile main script

backus_lsq: bgv6_leastsq.c libinterface.so libzcall.so
	$(CC) -o $@ $< $(CFLAGS) -lzcall -linterface $(ZLFLAGS) $(LFLAGS)

backus_spr: bgv6_spread.c libinterface.so libzcall.so
	$(CC) -o $@ $< $(CFLAGS) -lzcall -linterface $(ZLFLAGS) $(LFLAGS)

all:
	export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:lib/"
	@echo $LD_LIBRARY_PATH
