VERSION  = 0.0.1
BOOST_DIR = /home/koenig/boost-install
ROOT_DIR =$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
CC = nvcc
CFLAGS  = -isystem $(BOOST_DIR)/include -I $(ROOT_DIR)/include

CPP_TEST = $(ROOT_DIR)/tests
CPP_DIR  = $(ROOT_DIR)/src

CPP_SAMPLE1 = $(CPP_TEST)/sample1_unittest.cu 
CPP_DYN     = $(CPP_TEST)/dyn_unittest.cu
sample1: $(CPP_SAMPLE1) 
	$(CC) $(CFLAGS) -o sample1 $(CPP_SAMPLE1)
dyn: $(CPP_DYN)
	$(CC) $(CFLAGS) -o dyn $(CPP_DYN)

