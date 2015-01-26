#CXX=mpiCC
CXX=g++
#FLAGS=-O3 -Wall -Wno-unused-variable -Wno-unknown-pragmas  -DNDEBUG -DUSE_MPI
FLAGS=-O3 -Wall -Wno-unused-variable -Wno-unknown-pragmas  -DNDEBUG -static
LIBS=-lboost_serialization

include Makefile.common

.PHONY: depend
depend:
	gccmakedep -- $(CXXFLAGS) -- $(SRCS) 
# DO NOT DELETE
