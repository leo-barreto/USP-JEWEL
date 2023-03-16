all: jewel-usp

# path to LHAPDF library
LHAPDF_PATH := /sampa/archive/leonardo/Jets/lhapdf6/lib

FC := gfortran
FFLAGS := -O2

jewel-usp: jewel240mod-usp.o medium-hydro.o reader-hydro.o interpolate.o pythia6425mod-lhapdf6.o meix.o
	$(FC) -o $@ -L$(LHAPDF_PATH) $^ -lLHAPDF -lstdc++
clean:
	rm -f jewel*.o
	rm -f medium-*.o interpolate.o reader-*.o 
	rm -f pythia6425mod-lhapdf6.o meix.o
	rm -f *~

.PHONY: all
