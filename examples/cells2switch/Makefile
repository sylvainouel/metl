CXX=g++
#CXX=/usr/lib/gcc-snapshot/bin/g++
FLAGS= -Wall -Wno-unused-variable -Wno-unknown-pragmas -O3 -static 
LDFLAGS = -L../../lib -lm -lsprng -lgmp


include Makefile.common

.PHONY: depend
depend:
	gccmakedep -- $(CXXFLAGS) -- $(SRCS) 
# DO NOT DELETE
