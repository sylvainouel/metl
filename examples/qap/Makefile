#CXX=mpiCC
CXX=llvmg++
#FLAGS=-O3 -Wall -Wno-unused-variable -Wno-unknown-pragmas -ffast-math -static -DUSE_MPI -DUSE_PAR
FLAGS=-O3 -Wall -Wno-unused-variable -Wno-unknown-pragmas -ffast-math -static -fno-exceptions -DBOOST_NO_EXCEPTIONS 
LDFLAGS = -L../../lib -lm -lsprng -lgmp 
# -lboost_serialization

include Makefile.common

.PHONY: depend
depend:
	gccmakedep -- $(CXXFLAGS) -- $(SRCS) 
# DO NOT DELETE
