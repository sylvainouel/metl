.PHONY: all

all: c2s_seq_algo.cc c2s_mpi_algo.cc c2s_omp_algo.cc c2s_par_algo.cc qap_algo.cc

all_algo.cc: gen_algo.tcl
	tclsh gen_algo.tcl > $@

c2s_seq_algo.cc: all_algo.cc
	cat $< | grep -v qap_prob | grep -v mpi | grep -v omp > $@

c2s_mpi_algo.cc: all_algo.cc
	cat $< | grep -v qap_prob | grep mpi | grep -v omp > $@

c2s_omp_algo.cc: all_algo.cc
	cat $< | grep -v qap_prob | grep omp | grep -v mpi > $@

c2s_par_algo.cc: all_algo.cc
	cat $< | grep -v qap_prob | grep omp | grep mpi > $@

qap_algo.cc: all_algo.cc
	cat $< | grep -v cell2switch > $@

.PHONY: clean

clean:
	rm -f *.cc *~
