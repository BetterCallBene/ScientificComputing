VERSION  = 0.0.1
BOOST_DIR = /home/koenig/boost-install
ROOT_DIR =$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
CC = g++
CFLAGS  = -isystem $(BOOST_DIR)/include -I$(ROOT_DIR)/include

CPP_TEST = $(ROOT_DIR)/tests
CPP_DIR  = $(ROOT_DIR)/src

CPP_SAMPLE1 = $(CPP_TEST)/sample1_unittest.cpp 
CPP_DYN     = $(CPP_TEST)/dyn_unittest.cpp

sample1.o: $(CPP_SAMPLE1) 
	$(CC) $(CFLAGS) $(CPP_SAMPLE1) -o sample1 